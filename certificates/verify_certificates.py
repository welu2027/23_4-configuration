import json, glob, re, sys
from fractions import Fraction
term_re = re.compile(r'([+-]?)([^+-]+)')
def parse(txt, VI, NV):
    poly = {}
    for sign, body in term_re.findall(txt.replace(' ', '')):
        if not body: continue
        coef = Fraction(1); exps = [0] * NV
        for factor in body.split('*'):
            if '^' in factor:
                base, e = factor.split('^'); exps[VI[base]] += int(e)
            elif factor in VI: exps[VI[factor]] += 1
            else: coef *= Fraction(factor)
        if sign == '-': coef = -coef
        key = tuple(exps); poly[key] = poly.get(key, 0) + coef
    return {k: c for k, c in poly.items() if c != 0}
def mul(a, b):
    acc = {}
    for ka, ca in a.items():
        for kb, cb in b.items():
            k = tuple(x + y for x, y in zip(ka, kb)); acc[k] = acc.get(k, 0) + ca * cb
    return {k: v for k, v in acc.items() if v != 0}
def add(a, b):
    acc = dict(a)
    for k, v in b.items(): acc[k] = acc.get(k, 0) + v
    return {k: v for k, v in acc.items() if v != 0}
total = ok = 0
for fn in sorted(glob.glob('*_certificates.json')):
    d = json.load(open(fn)); vs = d['variables']; VI = {v: i for i, v in enumerate(vs)}; NV = len(vs); one = {tuple([0] * NV): Fraction(1)}
    for idx, s in d['structures'].items():
        lhs = {}
        for a, f in zip(s['certificate'], s['generators']): lhs = add(lhs, mul(parse(a, VI, NV), parse(f, VI, NV)))
        total += 1; ok += (lhs == one)
    print("%s: %d certificates checked, non-unit structures %s" % (fn, len(d['structures']), d['nonunit']))
d = json.load(open('cert_88_degenerate.json')); vs = json.load(open('23_B_4_4_3_3_certificates.json'))['variables']; VI = {v: i for i, v in enumerate(vs)}; NV = len(vs)
lhs = {}
for a, f in zip(d['coefficients'], d['generators']): lhs = add(lhs, mul(parse(a, VI, NV), parse(f, VI, NV)))
rhs = {tuple([0] * NV): Fraction(1)}
for _ in range(d['k']): rhs = mul(rhs, parse(d['forbidden_incidence_polynomial'], VI, NV))
print("structure 88: (%s)^%d in the ideal: %s" % (d['forbidden_incidence_polynomial'], d['k'], lhs == rhs))
print("Nullstellensatz certificates verified: %d of %d" % (ok, total))
sys.exit(0 if ok == total else 1)
