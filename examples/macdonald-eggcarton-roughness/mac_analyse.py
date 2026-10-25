#!/usr/bin/env python3
"""Mean statistics of the MacDonald et al. (JFM 2019) set-A runs (CaNS-IBM, constant pressure gradient, u_tau = 1).
usage: mac_analyse.py <run dir> <rough|smooth> <first step> [last step] [--out file.npz]
Averages the plane statistics written every iout1d steps in [first, last]; the two channel halves are folded.
Rough: intrinsic (fluid-only) means, z measured from the roughness mean height (k above the bottom crest-trough base).
"""
import sys, glob, re, numpy as np
run, case, s0 = sys.argv[1], sys.argv[2], int(sys.argv[3])
s1 = int(sys.argv[4]) if len(sys.argv) > 4 and not sys.argv[4].startswith('--') else 10**12
out = sys.argv[sys.argv.index('--out')+1] if '--out' in sys.argv else None
rough = case == 'rough'
retau, pr, k = 395., 0.7, 1./18.
g = np.array(open(f'{run}/data/grid.out').read().split(), float).reshape(-1, 5)  # 0, zf, zc, dzf, dzc, k = 0..n+1
zc, dzf = g[1:-1, 2], g[1:-1, 3]; n = zc.size; lz = g[-2, 1]
steps = sorted(int(re.search(r'_(\d+)\.out$', f).group(1)) for f in glob.glob(f'{run}/data/scal_01_mean_fld_*.out'))
steps = [s for s in steps if s0 <= s <= s1]
acc = None
for s in steps:
    S = np.loadtxt(f'{run}/data/scal_01_mean_fld_{s:07d}.out')[:, 1]
    if rough:
        t = np.loadtxt(f'{run}/data/turb_stats_ibm_{s:07d}.out')
        phu, phc, U, W = t[:, 1], t[:, 3], t[:, 4], t[:, 6]
        uu, vv, ww, uw = t[:, 7], t[:, 8], t[:, 9], t[:, 10]
    else:
        t = np.loadtxt(f'{run}/data/turb_stats_{s:07d}.out')
        phu = phc = np.ones(n); U, W = t[:, 1], t[:, 3]
        uu, vv, ww, uw = t[:, 4], t[:, 5], t[:, 6], t[:, 7]
    Th = np.where(phc > 0, S/np.where(phc > 0, phc, 1), 0.)
    # second moments about zero (so that the time average of the covariance includes the plane-mean fluctuation)
    rec = np.array([phu, phc, U, Th, uu+U**2, vv, ww, uw, phu*U, S])
    acc = rec if acc is None else acc + rec
nt = len(steps); a = acc/nt
phu, phc, U, Th, uu2, vv, ww, uw, Usup, Ssup = a
uu = uu2 - U**2
# fold the two halves (the grid is symmetric)
f = lambda q, sgn=1: 0.5*(q + sgn*q[::-1])
U, Th, uu, vv, ww, phu, phc = map(f, (U, Th, uu, vv, ww, phu, phc)); uw = f(uw, -1)
zoff = k if rough else 0.
zp = (zc - zoff)*retau                                  # from the mean height of the bottom wall
h = n//2
Q = np.sum(Usup*dzf)                                    # flow rate per unit span (superficial integral)
Ub = Q/2.                                               # fluid volume per unit plan area = 2h = 2
qw = Q/2.                                               # wall heat flux per wall = int(u*1) dz / 2 (source -u*dT/dx, dT/dx = -1)
Thtau = qw                                              # u_tau = 1
Thp = Th/Thtau
Thm = np.sum(Usup*Th*dzf)/Q/Thtau                      # mixed-mean temperature (+)
res = dict(nt=nt, steps=(steps[0], steps[-1]) if nt else None, Ub=Ub, Cf=2/Ub**2, Thm=Thm, Ch=1/(Ub*Thm),
           Ta=np.sum(phc*Th*dzf)/np.sum(phc*dzf)/Thtau, Reb=2*Ub*retau)
if rough:
    ik = np.argmin(abs(zc - 2*k))                       # crest
    res.update(Uk=np.interp(2*k, zc[:h], U[:h]), Thk=np.interp(2*k, zc[:h], Thp[:h]))
print(f"{run}: {case}, {nt} samples, steps {res['steps']}")
for kk in ['Ub', 'Cf', 'Thm', 'Ch', 'Ta', 'Reb'] + (['Uk', 'Thk'] if rough else []):
    print(f"  {kk:4s} = {res[kk]:.5g}")
if out:
    np.savez(out, zp=zp[:h], U=U[:h], Th=Thp[:h], uu=uu[:h], vv=vv[:h], ww=ww[:h], uw=uw[:h], phu=phu[:h], phc=phc[:h],
             **{kk: v for kk, v in res.items() if kk != 'steps'})
