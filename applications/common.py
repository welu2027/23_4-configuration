"""Shared helpers for the Section 8 ("Connections to other fields") computations.

Run every script from this folder (applications/). The two configurations are read
from the certificate files at the root of the repository.
"""
import json
import os
import random

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')


def load(name):
    """Incidence matrix (23 x 23, 0/1) of configuration 'A' or 'B'.

    For B the certificate lists points against their polars, so the matrix is
    symmetric with zero diagonal.
    """
    f = {'A': 'integer_23_4_certificate.json', 'B': 'geometric_23_4_certificate.json'}[name]
    with open(os.path.join(ROOT, f)) as fh:
        return np.array(json.load(fh)['incidence'], dtype=int)


def check_configuration(H):
    """True if H is the incidence matrix of a combinatorial (n_4) configuration."""
    O = H @ H.T
    np.fill_diagonal(O, 0)
    return set(H.sum(0)) == {4} and set(H.sum(1)) == {4} and O.max() <= 1


def switch_chain(H, rng, sweeps):
    """Random edge switches of the Levi graph.

    A switch replaces incidences (l1,p1),(l2,p2) by (l1,p2),(l2,p1). It is accepted only
    if the result is still a configuration (degrees 4, two lines meet at most once), so
    the degrees and the girth are preserved. Runs until sweeps * 92 switches are accepted.
    """
    H = H.copy()
    n = H.shape[0]
    E = [(i, j) for i in range(n) for j in range(n) if H[i, j]]
    acc = 0
    while acc < sweeps * len(E):
        a, b = rng.sample(range(len(E)), 2)
        (i1, j1), (i2, j2) = E[a], E[b]
        if i1 == i2 or j1 == j2 or H[i1, j2] or H[i2, j1]:
            continue
        H[i1, j1] = H[i2, j2] = 0
        H[i1, j2] = H[i2, j1] = 1
        ok = True
        for i in (i1, i2):
            ov = H @ H[i]
            ov[i] = 0
            if ov.max() > 1:
                ok = False
                break
        if ok:
            E[a] = (i1, j2)
            E[b] = (i2, j1)
            acc += 1
        else:
            H[i1, j2] = H[i2, j1] = 0
            H[i1, j1] = H[i2, j2] = 1
    return H


def sigma2(H):
    """Second largest singular value of the incidence matrix (expansion of the Levi graph)."""
    return float(np.linalg.svd(H.astype(float), compute_uv=False)[1])


def e_efficiency(H):
    """E-efficiency factor of the block design (points = treatments, lines = blocks).

    N^T N = 4I + M with M the collinearity adjacency, so the factor equals (16 - sigma2^2)/16.
    """
    return (16 - sigma2(H) ** 2) / 16


def collinearity(H):
    Pg = (H.T @ H > 0).astype(int)
    np.fill_diagonal(Pg, 0)
    return Pg


def min_relays(H):
    """Fewest relays over pairs of points (nodes) that are not collinear (share no key).

    A relay for a non-collinear pair is a point collinear with both of them.
    """
    Pg = collinearity(H)
    CN = Pg @ Pg
    n = len(Pg)
    return int(min(CN[i, j] for i in range(n) for j in range(i + 1, n) if not Pg[i, j]))
