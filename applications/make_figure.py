"""Two-panel figure comparing A and B with the random sample (made for Section 8).

Reads data/random_sample.json (from sample_random.py) and writes
connections_figure.pdf and connections_figure.png.
"""
import json
from collections import Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from common import load, sigma2, min_relays

with open('data/random_sample.json') as f:
    D = json.load(f)
s2 = np.array([d['sigma2'] for d in D])
rel = Counter(d['min_relays'] for d in D)
N = len(D)
A, B = load('A'), load('B')
sA, sB = sigma2(A), sigma2(B)
try:
    with open('data/anneal_best.json') as f:
        sBest = json.load(f)['sigma2']
except FileNotFoundError:
    sBest = None

cA, cB, grey, ink = '#1f77b4', '#d62728', '#a3aaa7', '#222222'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': '#555555', 'axes.linewidth': 0.6,
                     'xtick.color': '#555555', 'ytick.color': '#555555',
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.1, 3.1), gridspec_kw=dict(wspace=0.28))
for ax in (a1, a2):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

a1.hist(s2, bins=np.arange(2.54, 3.05, 0.01), color=grey, edgecolor='white', linewidth=0.6)
top = a1.get_ylim()[1]
for x, c, lab in [(sA, cA, f'A   {sA:.3f}'), (sB, cB, f'B   {sB:.3f}')]:
    a1.axvline(x, color=c, lw=1.4)
    a1.plot([x], [top * 0.93], 'o', ms=6, mfc=c, mec='black', mew=0.7, zorder=5)
    a1.text(x + 0.014, top * 0.93, lab, va='center', color=ink, fontsize=8.5,
            bbox=dict(fc='white', ec='none', pad=0.6))
if sBest is not None:
    a1.axvline(sBest, color='#555555', lw=0.9, ls=(0, (3, 2)))
    a1.text(sBest - 0.012, top * 0.55, f'best found\nby search\n{sBest:.3f}', ha='right',
            va='center', fontsize=7.5, color='#555555')
a1.set_xlim(2.38, 3.05)
a1.set_ylim(0, top)
a1.set_xlabel(r'second singular value $\sigma_2$ of the incidence matrix (smaller is better)', fontsize=8.5)
a1.set_ylabel('number of configurations', fontsize=8.5)
a1.set_title('(a) expansion of the Levi graph', fontsize=9.5)

ks = [3, 4, 5, 6]
pct = [100 * rel.get(k, 0) / N for k in ks]
a2.bar(ks, pct, width=0.62, color=grey, edgecolor='white', linewidth=0.6)
for k, p in zip(ks, pct):
    a2.text(k, p + 1.6, ('%.1f%%' % p) if p > 0 else '0%', ha='center', va='bottom', fontsize=8, color=ink)
for x, c, lab in [(min_relays(A), cA, 'A'), (min_relays(B), cB, 'B')]:
    y = pct[ks.index(x)] + 13
    a2.plot([x], [y], 'o', ms=7, mfc=c, mec='black', mew=0.7, zorder=5)
    a2.text(x + 0.17, y, lab, va='center', fontsize=9, color=ink)
a2.set_xticks(ks)
a2.set_ylim(0, 92)
a2.set_xlim(2.4, 6.6)
a2.set_xlabel('fewest relays over pairs of nodes with no common key', fontsize=8.5)
a2.set_ylabel('% of random configurations', fontsize=8.5)
a2.set_title('(b) key predistribution', fontsize=9.5)

fig.savefig('connections_figure.pdf', bbox_inches='tight')
fig.savefig('connections_figure.png', dpi=150, bbox_inches='tight')
print('wrote connections_figure.pdf and connections_figure.png')
