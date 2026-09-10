from itertools import product
def sums(sizes, target, maxcount=None):
    out = []
    def rec(i, rem, cur):
        if i == len(sizes):
            if rem == 0:
                out.append(tuple(cur))
            return
        s = sizes[i]
        for k in range(rem // s + 1):
            rec(i + 1, rem - k * s, cur + [k])
    rec(0, target, [])
    return out
def cyclic_case(m):
    if m % 2 == 1:
        special = m
    else:
        special = m // 2
    pt_sizes = [1, special, m]
    ln_sizes = [1, special, m]
    survivors = []
    for (a, b, c) in sums(pt_sizes, 23):
        if a > 1:
            continue
        for (a2, b2, c2) in sums(ln_sizes, 23):
            if a2 > 1:
                continue
            if a == 1:
                if b2 * special != 4 and not (special == 1):
                    if b2 * special != 4:
                        continue
                else:
                    if b2 * special != 4:
                        continue
            if a2 == 1:
                if b * special != 4:
                    continue
            stab = m // special
            if b > 0:
                ok = False
                for e1 in (0, 1):
                    for e2 in (0, 1):
                        if e1 == 1 and a2 == 0: continue
                        if e2 == 1 and b2 == 0: continue
                        if (4 - e1 - e2) % stab == 0 and (4 - e1 - e2) >= 0:
                            ok = True
                if not ok:
                    continue
            if b2 > 0:
                ok = False
                for e1 in (0, 1):
                    for e2 in (0, 1):
                        if e1 == 1 and a == 0: continue
                        if e2 == 1 and b == 0: continue
                        if (4 - e1 - e2) % stab == 0 and (4 - e1 - e2) >= 0:
                            ok = True
                if not ok:
                    continue
            if b * special > 4 and a2 == 1:
                continue
            survivors.append(((a, b, c), (a2, b2, c2)))
    return survivors
print("Cyclic groups C_m, m = 3..46 : orbit types surviving the counting step.")
print("  Expected: only single-orbit types with m = 23 or 46 (all 23 points on the axis line,")
print("  or one regular orbit of a rotation of order 23), which Proposition 6 excludes by the")
print("  line/conic argument, not by counting.")
for m in range(3, 47):
    s = cyclic_case(m)
    if s:
        assert m in (23, 46) and all(a == 0 and b + c == 1 for (a, b, c), _ in s), (m, s)
        print("  m =", m, "->", s, " [single orbit: excluded by the line/conic argument]")
print("  Counting step complete: every survivor is a single-orbit case with m in {23, 46}.")
groups = {
    'A4 (order 12)': [3, 4, 6, 12],
    'S4 (order 24)': [3, 4, 6, 12, 24],
    'A5 (order 60)': [6, 10, 15, 30, 60],
}
print("Polyhedral groups: orbit-size decompositions of 23 (points) and (lines) without local test:")
for name, sizes in groups.items():
    dec = sums(sizes, 23)
    print(" ", name, ":", len(dec), "decompositions", dec[:6], "..." if len(dec) > 6 else "")
print()
print("BUT A4, S4, A5 all contain C_3 (and S4 contains C_4); A5 contains C_5.")
print("Since no C_3, C_4, C_5 can act, none of these can.  Dihedral D_m (m>=3) contains C_m -> excluded.")
print("Remaining: C_2 and D_2 = C_2 x C_2 (Klein four-group).")
