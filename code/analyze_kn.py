import os
import json, glob, sys, itertools, math, os, numpy as np, pynauty, mpmath
KN_DIR = os.environ.get('KN_DIR', '/home/claude/kn')
os.makedirs(KN_DIR, exist_ok=True)
from fractions import Fraction as F
from classify_Vn import *
mpmath.mp.dps = 30
n = int(sys.argv[1])
def numeric(lay, sol):
    pts, lns = sym_coords(lay)
    P = np.array([[float(s) * (sol[var].real if var else 1.0) for (s, var) in pts[i]] for i in range(n)])
    L = np.array([[float(s) * (sol[var].real if var else 1.0) for (s, var) in lns[i]] for i in range(n)])
    return P, L
def levi(inc):
    g = pynauty.Graph(2 * n, directed=False, vertex_coloring=[set(range(n)), set(range(n, 2 * n))])
    g.set_adjacency_dict({p: [n + l for l in range(n) if inc[p][l]] for p in range(n)})
    gens, o1, o2, orbs, norb = pynauty.autgrp(g); return pynauty.certificate(g), int(round(o1 * 10 ** o2)), gens, pynauty.canon_label(g)
def group_elements(gens, N):
    els = {tuple(range(N))}; frontier = [tuple(range(N))]
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                x = tuple(h[g[i]] for i in range(N))
                if x not in els: els.add(x); new.append(x)
        frontier = new
    return els
def realized(A, B, perm):
    rows = []
    for k in range(n):
        p = A[k]; q = B[perm[k]]
        for (a, b) in ((0, 1), (1, 2), (0, 2)):
            row = np.zeros(9)
            for t in range(3): row[3 * a + t] += p[t] * q[b]; row[3 * b + t] -= p[t] * q[a]
            rows.append(row)
    sv = np.linalg.svd(np.array(rows), compute_uv=False); return sv[-1] < 1e-8 * sv[0]
sols = []; nstruct = 0; ncomplex = 0
for fn in sorted(glob.glob(KN_DIR + '/%d_*_results.json' % n)):
    data = json.load(open(fn.replace('_results.json', '.json'))); cs = data['case']; lay = Layout(n, cs[0], cs[1], cs[2], cs[3], cs[4])
    res = json.load(open(fn)); keys = sorted(data['structures']); nstruct += len(keys)
    for idx, key in enumerate(keys):
        r = res[key]
        if r.get('status') != 'ok': continue
        for s in r['sols']:
            if not s['real']: ncomplex += 1; continue
            sol = {k: complex(v[0], v[1]) for k, v in s['sol'].items()}
            P, L = numeric(lay, sol); inc = [[int(abs(P[i] @ L[j]) < 1e-9 * max(1, np.linalg.norm(P[i]) * np.linalg.norm(L[j]))) for j in range(n)] for i in range(n)]
            ok = all(sum(r_) == 4 for r_ in inc) and all(sum(inc[i][j] for i in range(n)) == 4 for j in range(n))
            cert, aut, gens, lab = levi(inc)
            sols.append({'case': fn.split('/')[-1].replace('_results.json', ''), 'idx': idx, 'vdim': r['vdim'], 'sol': sol, 'P': P, 'L': L, 'inc': inc, 'cert': cert, 'aut': aut, 'gens': gens, 'lab': lab, 'ok': ok})
print("n=%d: %d structures; real nondegenerate solutions: %d (all incidence-checked: %s); complex-only: %d" % (n, nstruct, len(sols), all(s['ok'] for s in sols), ncomplex))
classes = []
for s in sols:
    for cl in classes:
        t = cl[0]
        if t['cert'] != s['cert']: continue
        inv_t = {v: i for i, v in enumerate(t['lab'])}; base = [s['lab'][inv_t[v]] for v in range(n)]
        G = group_elements(t['gens'], 2 * n); perms = [[base[g[i]] for i in range(n)] for g in G if all(g[i] < n for i in range(n))]
        if any(realized(t['P'], s['P'], perm) for perm in perms): cl.append(s); break
    else: classes.append([s])
print("projective equivalence classes of real Klein-symmetric (%d_4)s: %d" % (n, len(classes)))
types = {}
for cl in classes: types.setdefault(cl[0]['cert'], []).append(cl)
print("combinatorial types:", len(types))
def minpoly(x):
    for deg in (1, 2, 4, 8):
        p = mpmath.findpoly(mpmath.mpf(x), deg, maxcoeff=10 ** 5)
        if p: return p
    return None
summary = []
for ti, (cert, cls) in enumerate(types.items()):
    s = cls[0][0]; G = group_elements(s['gens'], 2 * n); coll = sum(1 for g in G if all(g[i] < n for i in range(n)) and realized(s['P'], s['P'], list(g[:n])))
    incT = [[s['inc'][i][j] for i in range(n)] for j in range(n)]; selfdual = levi(incT)[0] == cert
    irr = [k for k in s['sol'] if abs(s['sol'][k].real - round(s['sol'][k].real)) > 1e-6]
    mps = sorted(set(str(minpoly(s['sol'][k].real)) for k in irr[:6]))
    print("type %d: %d projective class(es); |Aut(Levi)| = %d; collineations = %d; self-dual: %s; vdim %d; from %s#%d; minpolys of some coordinates: %s" % (ti, len(cls), s['aut'], coll, selfdual, s['vdim'], s['case'], s['idx'], mps))
    summary.append({'type': ti, 'classes': len(cls), 'aut': s['aut'], 'collineations': coll, 'selfdual': selfdual, 'vdim': s['vdim'], 'case': s['case'], 'idx': s['idx'], 'incidence': s['inc'], 'sol': {k: [v.real, v.imag] for k, v in s['sol'].items()}, 'minpolys': mps})
json.dump({'n': n, 'structures': nstruct, 'real_solutions': len(sols), 'complex_only': ncomplex, 'classes': len(classes), 'types': summary}, open(KN_DIR + '/summary_%d.json' % n, 'w'))
if n == 22:
    if os.path.exists('new_22_4_certificates.json'):
        known = {}
        for name, lst in json.load(open('new_22_4_certificates.json')).items(): known[levi(lst[0]['incidence'])[0]] = name
        for ti, cert in enumerate(types):
            if cert in known: print("  type %d = previously found %s" % (ti, known[cert]))
    else:
        print("  (no new_22_4_certificates.json; skipping comparison with the earlier run)")
if n == 23:
    for name, fn in (('A', 'integer_23_4_certificate.json'), ('B', 'geometric_23_4_certificate.json')):
        c = levi(json.load(open(fn))['incidence'])[0]
        print("  configuration %s found:" % name, c in types)
