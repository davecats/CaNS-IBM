#!/usr/bin/env python3
"""Compare the rough and smooth runs (npz from mac_analyse.py) with MacDonald et al. (JFM 2019), set A full span.
usage: mac_compare.py rough.npz smooth.npz [fig.png]"""
import sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
r, s = np.load(sys.argv[1]), np.load(sys.argv[2]); fig = sys.argv[3] if len(sys.argv) > 3 else 'mac_compare.png'
retau = 395.
# smooth profiles interpolated to the rough z+ (z+ > 0 above the mean height)
m = r['zp'] > 0
zp = r['zp'][m]
dU = np.interp(zp, s['zp'], s['U']) - r['U'][m]
dT = np.interp(zp, s['zp'], s['Th']) - r['Th'][m]
band = (zp > 50) & (zp < 0.9*retau)
c = np.argmin(abs(zp - retau))
print('                      present     MacDonald et al. (2019), table 1 set A (L_y+ = 1241)')
print(f"U_b+ rough        {float(r['Ub']):8.3f}     11.06 (Fig. 2a full-span profile; table U_bf+ = 10.5 is a composite-profile estimate)")
print(f"U_b+ smooth       {float(s['Ub']):8.3f}     17.7 (Fig. 2a full-span profile)")
print(f"dU+  centre       {dU[c]:8.3f}     6.7")
print(f"dU+  50<z+<355    {dU[band].mean():8.3f} (min {dU[band].min():.2f}, max {dU[band].max():.2f})")
print(f"dTh+ centre       {dT[c]:8.3f}     2.9")
print(f"dTh+ 50<z+<355    {dT[band].mean():8.3f} (min {dT[band].min():.2f}, max {dT[band].max():.2f})")
for kk in ['Cf', 'Ch', 'Thm', 'Reb']:
    print(f"{kk:4s} rough {float(r[kk]):10.5g}   smooth {float(s[kk]):10.5g}")
print(f"2Ch/Cf rough {2*float(r['Ch'])/float(r['Cf']):.4f}  smooth {2*float(s['Ch'])/float(s['Cf']):.4f}")
print(f"U_k+ {float(r['Uk']):.3f}  Th_k+ {float(r['Thk']):.3f}")
# digitised full-span curves of the paper's figure 2 (reference/ next to this script), if present
import os
pf = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reference', 'macdonald2019_fig2_digitised.dat')
paper = {}
if os.path.exists(pf):
    for line in open(pf):
        if line.startswith('#'): continue
        q, w, z, v = line.split(); paper.setdefault((q, w), []).append((float(z), float(v)))
f, ax = plt.subplots(1, 2, figsize=(11, 4.3))
for a, q, lab, d, ref in [(ax[0], 'U', '$U^+$', dU, 6.7), (ax[1], 'Th', r'$\Theta^+$', dT, 2.9)]:
    a.semilogx(s['zp'], s[q], 'k-', label='smooth')
    a.semilogx(r['zp'][m], r[q][m], 'k--', label='rough $k^+=22$')
    for w, mk in [('smooth', 'o'), ('rough', 's')]:
        pts = np.array(paper.get(('U' if q == 'U' else 'Theta', w), []))
        if pts.size: a.semilogx(pts[:, 0], pts[:, 1], mk, mfc='none', mec='C3', ms=4, mew=0.8, label=f'{w}, paper (fig. 2)')
    zz = np.logspace(1, np.log10(retau), 50)
    if q == 'U': a.semilogx(zz, np.log(zz)/0.4 + 5.0, ':', c='0.5', label=r'$\kappa_m^{-1}\ln z^+ + 5.0$')
    else: a.semilogx(zz, np.log(zz)/0.46 + 3.2, ':', c='0.5', label=r'$\kappa_h^{-1}\ln z^+ + 3.2$')
    a.axvline(21.9, c="0.8", lw=0.8)
    a.set_xlabel('$z^+$'); a.set_ylabel(lab); a.set_xlim(1, retau); a.set_ylim(0, None); a.legend(loc='upper left', fontsize=8)
    ins = a.inset_axes([0.62, 0.1, 0.35, 0.35])
    ins.semilogx(zp, d, 'k-'); ins.axhline(ref, c='C3', ls='--', lw=1); ins.set_xlim(1, retau)
    ins.set_title(r'smooth $-$ rough (dashed: paper, table 1)', fontsize=7); ins.tick_params(labelsize=7)
f.tight_layout(); f.savefig(fig, dpi=130)
print('wrote', fig)
