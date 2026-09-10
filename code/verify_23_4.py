from fractions import Fraction as F
import itertools, sys
def add(a, b): return tuple(x + y for x, y in zip(a, b))
def neg(a): return tuple(-x for x in a)
def mul(a, b):
    c = [F(0)] * 7
    for i in range(4):
        for j in range(4):
            c[i + j] += a[i] * b[j]
    for k in (6, 5, 4):
        if c[k]:
            c[k - 2] += c[k] / 2; c[k - 4] += c[k]; c[k] = F(0)
    return tuple(c[:4])
def const(q): return (F(q), F(0), F(0), F(0))
def is_zero(a): return all(x == 0 for x in a)
ONE, ZERO = const(1), const(0)
TAU = (F(0), F(1), F(0), F(0))
TAU2 = mul(TAU, TAU)
TAUINV = add(mul(TAU, TAU2), neg(mul(const(F(1, 2)), TAU)))
assert mul(TAU, TAUINV) == ONE, "1/tau formula"
assert add(mul(const(2), mul(TAU2, TAU2)), neg(add(TAU2, const(2)))) == ZERO, "minimal polynomial"
TAUINV2 = mul(TAUINV, TAUINV)
def fin(X, Z): return (X, ONE, Z)
points, labels = [(ONE, ZERO, ZERO)], ['O1']
for s in (1, -1):
    points.append(fin(ZERO, TAU if s == 1 else neg(TAU)));       labels.append('T(0,%+d tau)' % s)
for s in (1, -1):
    points.append(fin(ZERO, TAUINV if s == 1 else neg(TAUINV))); labels.append('T(0,%+d/tau)' % s)
for s in (1, -1):
    points.append(fin(const(2 * s), ZERO));                       labels.append('S(%+d,0)' % (2 * s))
orbits = [('C1', ONE, add(TAU, TAUINV)), ('C2', mul(const(2), TAU2), TAU),
          ('C3', ONE, add(TAU, neg(TAUINV))), ('C4', mul(const(2), TAUINV2), TAUINV)]
for name, X, Z in orbits:
    for sx in (1, -1):
        for sz in (1, -1):
            points.append(fin(X if sx == 1 else neg(X), Z if sz == 1 else neg(Z)))
            labels.append('%s(%s,%s)' % (name, '+' if sx == 1 else '-', '+' if sz == 1 else '-'))
assert len(points) == 23
def bform(P, Q):
    return add(neg(mul(P[0], Q[0])), mul(const(2), add(mul(P[1], Q[1]), mul(P[2], Q[2]))))
def main():
    n = 23
    for i, j in itertools.combinations(range(n), 2):
        P, Q = points[i], points[j]
        cross = [add(mul(P[1], Q[2]), neg(mul(P[2], Q[1]))), add(mul(P[2], Q[0]), neg(mul(P[0], Q[2]))),
                 add(mul(P[0], Q[1]), neg(mul(P[1], Q[0])))]
        if all(is_zero(c) for c in cross):
            print("FAIL: points %s and %s coincide" % (labels[i], labels[j])); sys.exit(1)
    inc = [[1 if is_zero(bform(points[i], points[j])) else 0 for j in range(n)] for i in range(n)]
    bad = [i for i in range(n) if sum(inc[i]) != 4]
    if bad:
        print("FAIL: line(s) with != 4 points:", [labels[i] for i in bad]); sys.exit(1)
    for i, j in itertools.combinations(range(n), 2):
        if sum(inc[i][k] and inc[j][k] for k in range(n)) > 1:
            print("FAIL: points %s, %s lie on two common lines" % (labels[i], labels[j])); sys.exit(1)
    absolute = sum(inc[i][i] for i in range(n))
    print("OK: 23 distinct points; each of the 23 polar lines contains exactly 4 of them;")
    print("    each point lies on exactly 4 lines; no two points share two lines; absolute points: %d." % absolute)
    print("    => a geometric (23_4) configuration, self-polar w.r.t. X^2 - 2Z^2 = 2, symmetric under (X,Z)->(+-X,+-Z).")
    print("\nIncidence table (line = polar of the named point : its four points)")
    for i in range(n):
        print("  polar of %-12s: %s" % (labels[i], ', '.join(labels[j] for j in range(n) if inc[i][j])))
    return 0
def verify_integral():
    from itertools import product
    pts = [(1, 0, 0)]
    pts += [(0, 1, s) for s in (1, -1)] + [(0, 3, 2 * s) for s in (1, -1)] + [(1, s, 0) for s in (1, -1)]
    for (a, b, c) in ((3, 2, 2), (2, 3, 2), (6, 1, 2), (6, 5, 2)):
        pts += [(a, e * b, d * c) for e, d in product((1, -1), repeat=2)]
    lns = [(1, 0, 0)]
    lns += [(0, 1, s) for s in (1, -1)] + [(0, 2, 3 * s) for s in (1, -1)] + [(1, 0, 3 * s) for s in (1, -1)]
    for (a, b, c) in ((1, 2, 2), (2, 6, 9), (2, 2, 1), (2, 2, 5)):
        lns += [(a, e * b, d * c) for e, d in product((1, -1), repeat=2)]
    assert len(pts) == 23 and len(lns) == 23
    def cross(u, v): return (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
    for i in range(23):
        for j in range(i + 1, 23):
            assert cross(pts[i], pts[j]) != (0, 0, 0), "coincident points"
            assert cross(lns[i], lns[j]) != (0, 0, 0), "coincident lines"
    inc = [[int(sum(a * b for a, b in zip(p, l)) == 0) for l in lns] for p in pts]
    assert all(sum(r) == 4 for r in inc), "a point is not on exactly 4 lines"
    assert all(sum(inc[p][l] for p in range(23)) == 4 for l in range(23)), "a line does not contain exactly 4 points"
    for i in range(23):
        for j in range(i + 1, 23):
            assert sum(inc[i][l] and inc[j][l] for l in range(23)) <= 1, "two points on two common lines"
    print("Part (A), integral configuration: 23 distinct points, 23 distinct lines, 4 points on every line,")
    print("4 lines through every point, no two points on two common lines -- geometric (23_4) verified with integer arithmetic.")
if __name__ == '__main__':
    rc = main()
    verify_integral()
    sys.exit(rc)
