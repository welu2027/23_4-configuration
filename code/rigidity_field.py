import json, itertools
from fractions import Fraction as F
def qmul(a, b):
    c = [F(0)] * 7
    for i in range(4):
        for j in range(4): c[i + j] += a[i] * b[j]
    for k in (6, 5, 4):
        c[k - 2] += c[k] / 2; c[k - 4] += c[k]; c[k] = F(0)
    return tuple(c[:4])
def qadd(a, b): return tuple(x + y for x, y in zip(a, b))
def qneg(a): return tuple(-x for x in a)
def qinv(a):
    M = [list(qmul(a, tuple(F(1) if i == k else F(0) for i in range(4)))) for k in range(4)]
    A = [[M[k][r] for k in range(4)] + [F(1) if r == 0 else F(0)] for r in range(4)]
    for c in range(4):
        piv = next(r for r in range(c, 4) if A[r][c] != 0); A[c], A[piv] = A[piv], A[c]
        A[c] = [x / A[c][c] for x in A[c]]
        for r in range(4):
            if r != c and A[r][c] != 0: A[r] = [x - A[r][c] * y for x, y in zip(A[r], A[c])]
    return tuple(A[r][4] for r in range(4))
def qzero(a): return all(x == 0 for x in a)
rcA = json.load(open('integer_23_4_certificate.json')); PA = [tuple(F(x) for x in p) for p in rcA['points']]; LA = [tuple(F(x) for x in l) for l in rcA['lines']]; incA = rcA['incidence']
rcB = json.load(open('geometric_23_4_certificate.json'))
PB = [tuple(tuple(F(c) for c in coord) for coord in pt) for pt in rcB['points']]
Bf = [-1, 2, 2]; LB = [tuple(tuple(F(Bf[k]) * c for c in pt[k]) for k in range(3)) for pt in PB]; incB = rcB['incidence']
def modp_rank(rows, p):
    M = [[x % p for x in r] for r in rows]; rank = 0; ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((r for r in range(rank, len(M)) if M[r][c]), None)
        if piv is None: continue
        M[rank], M[piv] = M[piv], M[rank]; inv = pow(M[rank][c], -1, p); M[rank] = [(x * inv) % p for x in M[rank]]
        for r in range(len(M)):
            if r != rank and M[r][c]: M[r] = [(x - M[r][c] * y) % p for x, y in zip(M[r], M[rank])]
        rank += 1
    return rank
def frac_mod(x, p): return (x.numerator * pow(x.denominator, -1, p)) % p
def jac_and_trivial(P, L, inc, tofield):
    J = []
    for i in range(23):
        for j in range(23):
            if inc[i][j]:
                row = [0] * 138
                for k in range(3): row[3 * i + k] = tofield(L[j][k]); row[69 + 3 * j + k] = tofield(P[i][k])
                J.append(row)
    assert len(J) == 92
    T = []
    for i in range(23):
        row = [0] * 138
        for k in range(3): row[3 * i + k] = tofield(P[i][k])
        T.append(row)
    for j in range(23):
        row = [0] * 138
        for k in range(3): row[69 + 3 * j + k] = tofield(L[j][k])
        T.append(row)
    for a in range(3):
        for b in range(3):
            row = [0] * 138
            for i in range(23): row[3 * i + a] = tofield(P[i][b])
            for j in range(23): row[69 + 3 * j + b] = -tofield(L[j][a])
            T.append(row)
    return J, T
p = 1000003
JA, TA = jac_and_trivial(PA, LA, incA, lambda x: frac_mod(x, p))
print("(A): rank of Jacobian mod p = %d ; rank of trivial motions (scalings + gl_3) mod p = %d ; 138 - 54 = 84" % (modp_rank(JA, p), modp_rank(TA, p)))
for q in (13, 17, 19, 47, 53, 59, 61, 67, 89, 101, 103, 137, 149, 151, 157, 179, 191, 193, 197):
    roots = [r for r in range(q) if (2 * r ** 4 - r ** 2 - 2) % q == 0]
    if roots: break
r = roots[0]; print("(B): reducing modulo p = %d with tau -> %d" % (q, r))
def evalB(coord): return sum(frac_mod(c, q) * pow(r, i, q) for i, c in enumerate(coord)) % q
JB, TB = jac_and_trivial(PB, LB, incB, evalB)
print("(B): rank of Jacobian mod p = %d ; rank of trivial motions mod p = %d" % (modp_rank(JB, q), modp_rank(TB, q)))
def det2(u, v, a, b): return qadd(qmul(u[a], v[b]), qneg(qmul(u[b], v[a])))
in17 = inQ = tot = 0; example = None
for j in range(23):
    pts = [PB[i] for i in range(23) if incB[i][j]]
    for (a, b) in ((0, 1), (0, 2), (1, 2)):
        if not qzero(det2(pts[0], pts[1], a, b)) and not qzero(det2(pts[0], pts[2], a, b)) and not qzero(det2(pts[1], pts[2], a, b)): break
    p1, p2, p3, p4 = pts
    num = qmul(det2(p1, p3, a, b), det2(p2, p4, a, b)); den = qmul(det2(p1, p4, a, b), det2(p2, p3, a, b))
    cr = qmul(num, qinv(den)); tot += 1
    if cr[1] == 0 and cr[3] == 0: in17 += 1
    if cr[1] == 0 and cr[2] == 0 and cr[3] == 0: inQ += 1
    if (cr[1] != 0 or cr[3] != 0) and example is None: example = (j, cr)
print("(B): cross-ratios of the 23 collinear quadruples: %d of %d lie in Q(sqrt17), %d in Q" % (in17, tot, inQ))
if example: print("     e.g. line %d: cross-ratio = %s + %s t + %s t^2 + %s t^3  (not in Q(sqrt17))" % (example[0], *example[1]))
in17 = inQ = tot = 0
for j in range(23):
    pts = [PA[i] for i in range(23) if incA[i][j]]
    d2 = lambda u, v, a, b: u[a] * v[b] - u[b] * v[a]
    for (a, b) in ((0, 1), (0, 2), (1, 2)):
        if d2(pts[0], pts[1], a, b) != 0 and d2(pts[0], pts[2], a, b) != 0 and d2(pts[1], pts[2], a, b) != 0: break
    p1, p2, p3, p4 = pts; cr = (d2(p1, p3, a, b) * d2(p2, p4, a, b)) / (d2(p1, p4, a, b) * d2(p2, p3, a, b)); tot += 1
print("(A): all %d cross-ratios rational (as they must be); sample values:" % tot, sorted(set(str(F(d2(pts[0], pts[2], a, b) * d2(pts[1], pts[3], a, b)) / F(d2(pts[0], pts[3], a, b) * d2(pts[1], pts[2], a, b))) for j in [0])))
