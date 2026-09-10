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

## Data
`integer_23_4_certificate.json`, `geometric_23_4_certificate.json`: coordinates and incidence matrices.
`exact_equivalence_119.json`: the projectivities of Proposition 4.4(b) over Q(sqrt17).

## Provenance
An abstract under this title was submitted to the 2027 Joint Mathematics Meetings (AMS Special Session
on Applied and Computational Algebraic Geometry, I) on September 7, 2026; the submission system's
confirmation for submission id 66522 is timestamped September 7, 2026, 8:59 PM. The abstract already
states both geometric (23_4) configurations, the resulting characterisation that geometric
(n_4) configurations exist if and only if n >= 18 and n != 19, the bound of four on the order of the
projective symmetry group of a geometric (23_4) configuration, and the exact Groebner-basis
classification of the Klein-symmetric case. Receipts are provided in `documentation/`.

An independent solution by W. Strinz was made public the following day, on September 8, 2026
(github.com/wstrinz/configuration-23-4). The two works were arrived at independently, and that
solution is acknowledged in a note in the paper.

## Documentation
`documentation/jmm-2027-abstract-submission-66522.pdf`: the abstract as recorded by the JMM submission
system.
`documentation/jmm-submission-screenshot-2026-09-09.png`: the confirmation e-mail of September 7, 2026.
