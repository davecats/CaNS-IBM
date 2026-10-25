#!/usr/bin/env python3
"""Local wall shear stress and wall pressure of the rough channel, from the per-cell wall force (FV-fric_*).
usage: mac_wall.py <run dir> [--nx 512 --ny 256 --nz 480 --nlx 16 --nly 8] [--out wall.npz] [--dtlog <log>]
Run with ibm_iout_wall > 0. For every snapshot, the force of the fluid on the wall in each cut cell is split into its
normal (pressure) and tangential (viscous) parts. They are summed per (x, y) column, over the 128 roughness elements
(phase average) and over both walls: the upper wall is the lower one shifted by half a wavelength in x, mirrored in z.
Stresses in wall units (u_tau = 1, rho = 1); the wall pressure is p_w - <p>, from the normal stress F.n/A.
Outputs (npz): element maps on the (x, y) columns of one element, area-weighted over the cut cells of each column:
  tx, ty, tz  tangential (viscous) stress vector;  tm  its magnitude;  pw  wall pressure
  fx_p, fx_v  streamwise pressure and viscous force per unit plan area of the column;  zw  wall height above the mean
plus the instantaneous maps of the lower wall at the last snapshot (inst_tx, inst_pw) and the global totals.
--dtlog <log of the run>: for output of builds before e7e9133, where FV-fric used the dt of the next step on
steps that are multiples of icheck; the forces are rescaled by r = dt_next/dt. This leaves an error
(r - 1)*Pi*V_solid per cell (V_solid: solid volume of the cut cell, Pi the driving force), about 1e-3 of the drag
for r - 1 = 0.1-0.2 in the reference run. The lower-wall total is printed next to the exact wall-stress map
(ibm_wall_stress_<step>.out) when it exists.
"""
import sys, glob, re, argparse, numpy as np
ap = argparse.ArgumentParser()
ap.add_argument('run'); ap.add_argument('--nx', type=int, default=512); ap.add_argument('--ny', type=int, default=256)
ap.add_argument('--nz', type=int, default=480); ap.add_argument('--nlx', type=int, default=16); ap.add_argument('--nly', type=int, default=8)
ap.add_argument('--first', type=int, default=0); ap.add_argument('--out', default='wall.npz'); ap.add_argument('--dtlog')
o = ap.parse_args()
load = lambda f: np.fromstring(open(f).read(), sep=' ').reshape(-1, 26)    # fast reader of the text files
k, lam = 1./18., 2*np.pi/16.                                    # lambda = 7.07k
mx, my = o.nx//o.nlx, o.ny//o.nly                                # columns per element
Lx, Ly = o.nlx*lam, o.nly*lam
dx, dy = Lx/o.nx, Ly/o.ny
steps = sorted({int(re.search(r'_s_(\d+)\.out$', f).group(1)) for f in glob.glob(f'{o.run}/data/FV-fric_*_s_*.out')})
steps = [s for s in steps if s >= o.first]
dts = {}                                                         # dt set at each check, by step
if o.dtlog:
    st = 0
    for line in open(o.dtlog):
        if 'Time step #' in line: st = int(line.split()[3])
        elif 'dt_cfl =' in line: dts[st] = float(line.split()[-1])
def scale(s):
    if not o.dtlog or s not in dts: return 1.
    return dts[s]/dts[max(t for t in dts if t < s)]
acc = {q: np.zeros((mx, my)) for q in ['A', 'ftx', 'fty', 'ftz', 'ftm', 'fn', 'fpx', 'fvx', 'n']}
tot = np.zeros(3); nsnap = 0
for s in steps:
    a = np.concatenate([load(f) for f in glob.glob(f'{o.run}/data/FV-fric_*_s_{s:07d}.out')])
    r = scale(s); a[:, 7:19] *= r; a[:, 19:26] *= r                # force columns (not the area and normal)
    i, j, kk = a[:, 0].astype(int) - 1, a[:, 1].astype(int) - 1, a[:, 2].astype(int)
    up = kk > o.nz//2
    ie = np.where(up, (i + mx//2) % mx, i % mx); je = j % my
    for q, v in [('A', a[:, 3]), ('ftx', a[:, 15]), ('fty', a[:, 16]), ('ftz', np.where(up, -a[:, 17], a[:, 17])),
                 ('ftm', a[:, 14]), ('fn', a[:, 10]), ('fpx', a[:, 11]), ('fvx', a[:, 15]), ('n', np.ones(len(a)))]:
        np.add.at(acc[q], (ie, je), v)
    tot += [a[:, 7].sum(), a[:, 11].sum(), a[:, 15].sum()]
    nsnap += 1
    lo = ~up                                                     # instantaneous maps of the lower wall
    Ai = np.zeros((o.nx, o.ny)); txi = np.zeros((o.nx, o.ny)); fni = np.zeros((o.nx, o.ny))
    np.add.at(Ai, (i[lo], j[lo]), a[lo, 3]); np.add.at(txi, (i[lo], j[lo]), a[lo, 15]); np.add.at(fni, (i[lo], j[lo]), a[lo, 10])
    m = f'{o.run}/data/ibm_wall_stress_{s:07d}.out'
    chk = f', map {np.loadtxt(m)[:, 2].mean():.5f}' if glob.glob(m) else ''
    print(f'step {s}: {len(a)} cells, dt factor {r:.4f}, lower-wall drag per plan area {a[~up, 7].sum()/(Lx*Ly):.5f}{chk}; '
          f'per wall {a[:, 7].sum()/(2*Lx*Ly):.5f} (pressure {a[:, 11].sum()/(2*Lx*Ly):.5f})', flush=True)
nw = 2*o.nlx*o.nly*nsnap                                         # samples per element column
A = acc['A']
res = dict(tx=acc['ftx']/A, ty=acc['fty']/A, tz=acc['ftz']/A, tm=acc['ftm']/A, pw=-acc['fn']/A,
           fx_p=acc['fpx']/nw/(dx*dy), fx_v=acc['fvx']/nw/(dx*dy), area=A/nw/(dx*dy),
           x=(np.arange(mx) + 0.5)*dx/lam, y=(np.arange(my) + 0.5)*dy/lam,
           inst_tx=txi/np.where(Ai > 0, Ai, 1), inst_pw=-fni/np.where(Ai > 0, Ai, 1), inst_step=steps[-1] if steps else -1,
           nsnap=nsnap, steps=np.array(steps), drag=tot[0]/nsnap/(2*Lx*Ly), drag_p=tot[1]/nsnap/(2*Lx*Ly),
           drag_v=tot[2]/nsnap/(2*Lx*Ly))
xx, yy = np.meshgrid((np.arange(mx) + 0.5)*dx, (np.arange(my) + 0.5)*dy, indexing='ij')
res['zw'] = np.cos(2*np.pi*xx/lam)*np.cos(2*np.pi*yy/lam)        # wall height above the mean, in units of k
np.savez(o.out, **res)
print(f'{nsnap} snapshots; drag per wall and plan area {res["drag"]:.5f}: pressure {res["drag_p"]:.5f} '
      f'({res["drag_p"]/res["drag"]:.3f}), viscous {res["drag_v"]:.5f}; wrote {o.out}')
