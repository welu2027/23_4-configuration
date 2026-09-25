# Closing the existence problem for geometric (n_4) configurations

Two geometric (23_4) configurations, the symmetry theorem, and the classification of the Klein-symmetric case.

Abstract submitted September 7, 2026; note added September 9, 2026.

## Verify the main theorem
    python3 code/check_23_4.py
    python3 code/verify_23_4.py
    lean lean/Config23_4.lean
    lean lean/Config23_4_B.lean

## Verify the classification
    cd certificates && python3 verify_certificates.py
checks the 4018 Nullstellensatz certificates (one per empty polynomial system) and the certificate for the
degenerate structure, with exact rational arithmetic and no computer-algebra system.

## Reproduce the classification
    python3 code/classify_V_general.py 4 3 3 0 75
    python3 code/solve_general.py <case>.json 0 2000
    python3 code/dedup_audit.py <case>
    python3 code/exact_equiv.py
    python3 code/rigidity_field.py
Requirements for reproduction: Python 3.12 with pynauty, numpy, mpmath, sympy; Singular; PARI/GP.

## Connections to other fields (Section 8)
    cd applications
    python3 exact_B.py
    python3 sample_random.py
    python3 anneal_search.py
    python3 make_figure.py
proves the expansion and relay properties of (PB, LB) exactly, and reproduces the comparison with 20,000
random combinatorial (23_4) configurations. Requirements: numpy, networkx, sympy, matplotlib.
See `applications/README.md`.

## Data
`integer_23_4_certificate.json`, `geometric_23_4_certificate.json`: coordinates and incidence matrices.
`exact_equivalence_119.json`: the projectivities of Proposition 4.4(b) over Q(sqrt17).

## Provenance
After a preprint of this work was posted on Zenodo and submitted to the
Joint Mathematics Meetings 2027, a (23_4) configuration, which coincides with (PA, LA) up to the
coordinate interchange y ↔ z, was found independently by W. Strinz, and is accurately cited in the paper.
