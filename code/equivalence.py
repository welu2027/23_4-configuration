import json, cmath, itertools
import numpy as np, pynauty
from classify_polarV import *
def real_points(sol):
    thu = cmath.phase(sol['u1']) % cmath.pi
    thv = cmath.phase(sol['v1']) % cmath.pi
    mu = cmath.exp(-1j * thu); nu = cmath.exp(-1j * thv)
    pts = sym_points()
    env = dict(sol)
    out = []
    for i in range(23):
        P = [complex(eval(c.replace('^', '**'), {}, env)) for c in pts[i]]
        P = [P[0], mu * P[1], nu * P[2]]
        piv = [z for z in P if abs(z) > 1e-12][0]
        P = [z / piv for z in P]
        assert max(abs(z.imag) for z in P) < 1e-9, P
        v = np.array([z.real for z in P])
        out.append(v / np.linalg.norm(v))
    return np.array(out)
def levi_canon_labeling(A):
    n = 23
    g = pynauty.Graph(2 * n, directed=False, vertex_coloring=[set(range(n)), set(range(n, 2 * n))])
    g.set_adjacency_dict({i: [n + k for k in range(n) if A[i][k]] for i in range(n)})
    lab = pynauty.canon_label(g)
    gens, o1, o2, orbs, norb = pynauty.autgrp(g)
    return lab, gens, o1 * 10 ** o2
def group_elements(gens, n):
    ident = tuple(range(n))
    elems = {ident}
    frontier = [ident]
    gens = [tuple(g) for g in gens]
    while frontier:
        new = []
        for e in frontier:
            for g in gens:
                h = tuple(g[e[i]] for i in range(n))
                if h not in elems:
                    elems.add(h); new.append(h)
        frontier = new
    return elems
def proj_map(src, dst):
    A = src[:4].T; B = dst[:4].T
    try:
        la = np.linalg.solve(A[:, :3], A[:, 3]); lb = np.linalg.solve(B[:, :3], B[:, 3])
    except np.linalg.LinAlgError:
        return None
    if min(abs(la)) < 1e-9 or min(abs(lb)) < 1e-9: return None
    Ma = A[:, :3] * la; Mb = B[:, :3] * lb
    return Mb @ np.linalg.inv(Ma)
def equivalent(R1, R2, perms, n=23):
    for perm in perms:
        dst = R2[[perm[i] for i in range(n)]]
        for quad in itertools.combinations(range(n), 4):
            sub = R1[list(quad)]
            if abs(np.linalg.det(sub[:3])) < 1e-6 or abs(np.linalg.det(sub[1:])) < 1e-6 or\
               abs(np.linalg.det(sub[[0, 1, 3]])) < 1e-6 or abs(np.linalg.det(sub[[0, 2, 3]])) < 1e-6:
                continue
            M = proj_map(sub, dst[list(quad)])
            if M is None: continue
            ok = True
            for i in range(n):
                w = M @ R1[i]; w = w / np.linalg.norm(w)
                t = dst[i]
                if min(np.linalg.norm(w - t), np.linalg.norm(w + t)) > 1e-8:
                    ok = False; break
            if ok: return True
            break
    return False
if __name__ == '__main__':
    structs = json.load(open('structures_diag.json'))
    res = json.load(open('results_diag.json'))
    reals = []
    for tag, r in res.items():
        if r['status'] != 'ok' or r['nondeg'] == 0: continue
        idx = int(tag.split('_')[1]); A = structs[idx]['A']
        for s in r['sols']:
            sol = {k: complex(v[0], v[1]) for k, v in s['sol'].items()}
            reals.append({'tag': tag, 'A': A, 'sol': sol, 'R': real_points(sol)})
    print("nondegenerate solutions:", len(reals))
    lab0, gens0, order0 = levi_canon_labeling(reals[0]['A'])
    print("automorphism group order of the Levi graph:", order0)
    G = group_elements(gens0, 46)
    print("group elements generated:", len(G))
    inv0 = {v: i for i, v in enumerate(lab0)}
    classes = []
    for r in reals:
        lab, _, _ = levi_canon_labeling(r['A'])
        base = [lab[inv0[v]] for v in range(23)]
        assert all(b < 23 for b in base)
        placed = False
        for cl in classes:
            r0 = cl[0]
            lab_c, _, _ = levi_canon_labeling(r0['A']); inv_c = {v: i for i, v in enumerate(lab_c)}
            base_c = [lab[inv_c[v]] for v in range(23)]
            perms = []
            for g in G:
                pass
            _, gens_c, _ = levi_canon_labeling(r0['A'])
            Gc = group_elements(gens_c, 46)
            perms = [[base_c[g[i]] for i in range(23)] for g in Gc if all(g[i] < 23 for i in range(23))]
            if equivalent(r0['R'], r['R'], perms):
                cl.append(r); placed = True; break
        if not placed:
            classes.append([r])
    print("projective equivalence classes among real realizations:", len(classes))
    for k, cl in enumerate(classes):
        r = cl[0]
        print("class %d: size %d, e.g. %s  t1=%s t2=%s" % (k, len(cl), r['tag'], r['sol']['t1'], r['sol']['t2']))
    json.dump([{'tag': r['tag'], 'sol': {k: [v.real, v.imag] for k, v in r['sol'].items()}} for cl in classes for r in cl[:1]],
              open('real_classes.json', 'w'), indent=1)
