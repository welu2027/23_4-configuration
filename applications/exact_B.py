"""Exact proofs of the two properties of the self-polar configuration B (Section 8).

1. Expansion. B's incidence matrix M is symmetric (points listed against their polars),
   so its singular values are the absolute values of its eigenvalues. We factor the
   characteristic polynomial of M over Q. sigma_2 is the largest root of
   x^6 - 11x^4 + 31x^2 - 13, i.e. sigma_2^2 is the largest root of y^3 - 11y^2 + 31y - 13.

2. Relays. The automorphism group of the Levi graph of B (points and lines kept apart)
   has order 8. It splits the 115 pairs of non-collinear points into 18 orbits, and we
   list one representative per orbit with its number of relays (6, 7 or 8). Since relay
   counts are invariant under automorphisms, every non-collinear pair has at least 6.
"""
from collections import Counter

import networkx as nx
import numpy as np
import sympy as sp
from networkx.algorithms.isomorphism import GraphMatcher

from common import load, collinearity


def levi(H):
    G = nx.Graph()
    n = H.shape[0]
    G.add_nodes_from([('l', i) for i in range(n)] + [('p', j) for j in range(n)])
    G.add_edges_from((('l', i), ('p', j)) for i in range(n) for j in range(n) if H[i, j])
    return G


if __name__ == '__main__':
    M = load('B')
    assert (M == M.T).all() and np.trace(M) == 0

    # 1. exact spectrum
    x = sp.symbols('x')
    cp = sp.factor_list(sp.Matrix(M.tolist()).charpoly(x).as_expr())
    print('Characteristic polynomial of B (eigenvalues; singular values are their absolute values):')
    for fac, mult in cp[1]:
        roots = sorted(float(sp.re(r)) for r in sp.Poly(fac, x).nroots())
        print(f'  ({sp.expand(fac)})^{mult}   roots {[round(r, 4) for r in roots]}')
    y = sp.symbols('y')
    cubic = y ** 3 - 11 * y ** 2 + 31 * y - 13
    s2sq = max(float(r) for r in sp.Poly(cubic, y).nroots() if abs(sp.im(r)) < 1e-12)
    sv = np.sort(np.linalg.svd(M.astype(float), compute_uv=False))[::-1]
    print(f'largest root of y^3-11y^2+31y-13 = {s2sq:.6f}, sqrt = {s2sq ** 0.5:.6f}, '
          f'numerical sigma_2 = {sv[1]:.6f}')
    assert abs(s2sq ** 0.5 - sv[1]) < 1e-9

    # 2. relay orbits
    G = levi(M)
    auts = [a for a in GraphMatcher(G, G).isomorphisms_iter() if all(a[v][0] == v[0] for v in G)]
    perms = [tuple(a[('p', j)][1] for j in range(23)) for a in auts]
    print(f'\nautomorphism group of the Levi graph (sides kept apart): order {len(perms)}')
    Pg = collinearity(M)
    CN = Pg @ Pg
    pairs = [(i, j) for i in range(23) for j in range(i + 1, 23) if not Pg[i, j]]
    seen, orbits = set(), []
    for (i, j) in pairs:
        if (i, j) in seen:
            continue
        orb = {tuple(sorted((g[i], g[j]))) for g in perms}
        seen |= orb
        orbits.append(((i, j), len(orb), int(CN[i, j])))
    print(f'non-collinear pairs: {len(pairs)}, orbits: {len(orbits)}')
    print('representative   orbit size   relays')
    for rep, size, r in orbits:
        print(f'  points {rep[0]:2d},{rep[1]:2d}        {size}          {r}')
    assert sum(s for _, s, _ in orbits) == len(pairs)
    print('relay counts over orbits:', dict(sorted(Counter(r for *_, r in orbits).items())))
    print('minimum relay count of B:', min(r for *_, r in orbits))
