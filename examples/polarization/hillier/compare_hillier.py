#!/usr/bin/env python
'''
Compare the Sirocco runs with Hillier (1994) Figs 2-5; writes
compare_hillier.txt and compare_hillier.png.

usage: compare_hillier.py
'''
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import hillier_pol as hp
import hillier_ref as ref

lines = []
fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))

# P vs inclination for each chi0
cmap = plt.get_cmap('viridis')
for k, (chi0, href) in enumerate(ref.P_INC.items()):
    root = f'hx{chi0}'
    b = hp.band(root + '.spec')
    q = np.array([100 * Q / I for o, I, Q, U in b])
    u = np.array([100 * U / I for o, I, Q, U in b])
    c = cmap(k / (len(ref.P_INC) - 1))
    ax[0].plot(ref.INC, q, 'o', color=c, ms=5, label=f'{chi0}')
    hx = [x for x, y in zip(ref.INC, href) if y is not None]
    hy = [y for y in href if y is not None]
    ax[0].plot(hx, hy, '-', color=c, lw=1)
    lines.append(f'chi0 = {chi0}  (tau_av = {2.888*chi0:.3f})')
    lines.append('   i      Sirocco q    u      Hillier   diff')
    for i, qq, uu, hh in zip(ref.INC, q, u, href):
        hs = f'{hh:+6.2f}  {qq-hh:+6.2f}' if hh is not None else '   --'
        lines.append(f'  {i:6.2f}  {qq:+7.3f}  {uu:+7.3f}   {hs}')
    # q vs chi0 at each inclination (Fig 2 form) collected below
ax[0].set_xlabel('inclination (deg)')
ax[0].set_ylabel('q (%)  (> 0: parallel to axis)')
ax[0].set_title('pure scattering: points Sirocco, lines Hillier Fig 3')
ax[0].legend(title='chi0', fontsize=8, ncol=2)
ax[0].axhline(0, color='0.7', lw=0.5)

# Fig 2 form: q vs chi0 at 90 and 45 deg
chis = list(ref.P_INC)
for i_idx, mk in [(7, 'o'), (4, 's'), (3, '^'), (1, 'v')]:
    qs = [100 * hp.band(f'hx{c}.spec')[i_idx][2] / hp.band(f'hx{c}.spec')[i_idx][1] for c in chis]
    hs = [ref.P_INC[c][i_idx] for c in chis]
    l, = ax[1].plot(chis, qs, mk, label=f'{ref.INC[i_idx]}')
    ax[1].plot(chis, hs, '-', color=l.get_color(), lw=1)
cc = np.geomspace(0.01, 0.06, 20)
ax[1].plot(cc, -66.6 * cc, 'k:', lw=1, label='thin, 90')
ax[1].set_xscale('log')
ax[1].set_ylim(-4, 2)
ax[1].set_xlabel('chi0  (tau_av = 2.888 chi0)')
ax[1].set_title('Hillier Fig 2 form (lines Hillier)')
ax[1].legend(title='i', fontsize=8)
ax[1].axhline(0, color='0.7', lw=0.5)

# albedo
for root, ang, mk in [('hal90', 90.0, 'o'), ('hal56', 56.25, 's')]:
    res = hp.albedo(root)
    a = np.array([r[0] for r in res])
    q = np.array([r[1] for r in res])
    l, = ax[2].plot(a, q, mk, ms=5, label=f'Sirocco {ang}')
    ax[2].plot(ref.ALB, ref.P_ALB[ang], '-', color=l.get_color(), lw=1, label=f'Hillier {ang}')
    lines.append(f'albedo, i = {ang}, chi0 = 1')
    lines.append('   a     Sirocco q    u     Hillier   diff')
    for aa, h in zip(ref.ALB, ref.P_ALB[ang]):
        j = int(np.argmin(np.abs(a - aa)))
        lines.append(f'  {aa:4.2f}  {res[j][1]:+7.3f}  {res[j][2]:+7.3f}   {h:+6.2f}  {res[j][1]-h:+6.2f}')
ax[2].set_xlabel('albedo')
ax[2].set_title('chi0 = 1: Hillier Figs 4 (90) and 5 (56.25)')
ax[2].legend(fontsize=8)
ax[2].axhline(0, color='0.7', lw=0.5)

fig.tight_layout()
fig.savefig('compare_hillier.png', dpi=110)
open('compare_hillier.txt', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
