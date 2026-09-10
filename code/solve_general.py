import json, sys, subprocess, re, time, itertools, cmath, glob
from classify_V_general import *
def sym_coords(lay):
    pts = {0: ((1, None), (0, None), (0, None))}
    pts[lay.T(1, 0)] = ((0, None), (1, None), (1, None)); pts[lay.T(1, 1)] = ((0, None), (1, None), (-1, None))
    pts[lay.T(2, 0)] = ((0, None), (1, None), (1, 't2')); pts[lay.T(2, 1)] = ((0, None), (1, None), (-1, 't2'))
    if lay.b_axis == 3:
        pts[lay.B(0)] = ((1, None), (1, None), (0, None)); pts[lay.B(1)] = ((1, None), (-1, None), (0, None))
    else:
        pts[lay.B(0)] = ((1, None), (0, None), (1, None)); pts[lay.B(1)] = ((1, None), (0, None), (-1, None))
    for j in range(1, lay.c + 1):
        for al in V:
            pts[lay.P(j, al)] = ((1, None), (al[0] * al[1], 'u%d' % j), (al[0] * al[2], 'v%d' % j))
    lns = {0: ((1, None), (0, None), (0, None))}
    lns[lay.Tp(1, 0)] = ((0, None), (1, None), (-1, None)); lns[lay.Tp(1, 1)] = ((0, None), (1, None), (1, None))
    lns[lay.Tp(2, 0)] = ((0, None), (1, 't2'), (-1, None)); lns[lay.Tp(2, 1)] = ((0, None), (1, 't2'), (1, None))
    for k, cen in enumerate(lay.bp):
        var = 'S%d' % k
        if cen == 3:
            lns[lay.Bp(k, 0)] = ((1, None), (1, var), (0, None)); lns[lay.Bp(k, 1)] = ((1, None), (-1, var), (0, None))
        else:
            lns[lay.Bp(k, 0)] = ((1, None), (0, None), (1, var)); lns[lay.Bp(k, 1)] = ((1, None), (0, None), (-1, var))
    for a in range(1, lay.cp + 1):
        for al in V:
            lns[lay.G(a, al)] = ((1, None), (al[0] * al[1], 'm%d' % a), (al[0] * al[2], 'n%d' % a))
    return pts, lns
def dotpoly(P, L):
    terms = {}
    for (s1, x1), (s2, x2) in zip(P, L):
        if s1 == 0 or s2 == 0: continue
        mono = tuple(sorted(v for v in (x1, x2) if v is not None))
        terms[mono] = terms.get(mono, 0) + s1 * s2
    terms = {m: c for m, c in terms.items() if c != 0}
    if not terms: return None
    items = sorted(terms.items())
    if items[0][1] < 0: items = [(m, -c) for m, c in items]
    return tuple(items)
def poly_str(p):
    out = []
    for mono, c in p:
        s = ('%d' % c) if not mono else (('' if c == 1 else '-' if c == -1 else '%d*' % c) + '*'.join(mono))
        out.append(s)
    return '+'.join(out).replace('+-', '-')
def variables(lay):
    vs = ['t2'] + ['u%d' % j for j in range(1, lay.c + 1)] + ['v%d' % j for j in range(1, lay.c + 1)]
    vs += ['S%d' % k for k in range(len(lay.bp))]
    vs += ['m%d' % a for a in range(1, lay.cp + 1)] + ['n%d' % a for a in range(1, lay.cp + 1)]
    return vs
def build(lay, pl):
    pts, lns = sym_coords(lay)
    eqs, ineqs = set(), set()
    for p in range(23):
        for l in range(23):
            poly = dotpoly(pts[p], lns[l])
            if pl[p] >> l & 1:
                if poly is None: continue
                if len(poly) == 1 and poly[0][0] == (): return None, None
                eqs.add(poly)
            else:
                if poly is None: return None, None
                if len(poly) == 1 and poly[0][0] == (): continue
                ineqs.add(poly)
    return sorted(eqs), sorted(ineqs)
def singular_script(lay, eqs, tag):
    vs = variables(lay)
    s = 'ring r = 0,(%s),dp;\n' % ','.join(vs)
    s += 'ideal I = %s;\nideal J = std(I);\nint d0 = dim(J);\nif (d0 >= 0) {\n' % ','.join(poly_str(e) for e in eqs)
    for v in vs: s += '  J = sat(J, ideal(%s));\n' % v
    s += '  J = sat(J, ideal(t2^2-1));\n'
    for a, b in itertools.combinations(range(1, lay.c + 1), 2):
        s += '  J = sat(J, ideal(u%d^2-u%d^2, v%d^2-v%d^2));\n' % (a, b, a, b)
    for a, b in itertools.combinations(range(1, lay.cp + 1), 2):
        s += '  J = sat(J, ideal(m%d^2-m%d^2, n%d^2-n%d^2));\n' % (a, b, a, b)
    for k1, k2 in itertools.combinations(range(len(lay.bp)), 2):
        if lay.bp[k1] == lay.bp[k2]: s += '  J = sat(J, ideal(S%d^2-S%d^2));\n' % (k1, k2)
    s += '  J = std(J);\n}\nint d = dim(J);\n'
    s += 'print("RES %s dim0=" + string(d0) + " dim=" + string(d) + " vdim=" + string(vdim(J)));\n' % tag
    s += 'if (d == 0 and vdim(J) > 0) { def R = solve(J, 40, 0, "nodisplay"); setring R; print("SOLS %s"); print(SOL); print("ENDSOLS"); kill R; }\nkill r;\n' % tag
    return s
def parse_output(out, vs):
    res = {}
    for m in re.finditer(r'RES (\S+) dim0=(-?\d+) dim=(-?\d+) vdim=(-?\d+)', out):
        res[m.group(1)] = {'dim0': int(m.group(2)), 'dim': int(m.group(3)), 'vdim': int(m.group(4)), 'sols': []}
    for m in re.finditer(r'SOLS (\S+)\n(.*?)ENDSOLS', out, flags=re.S):
        vals = [l.strip() for l in m.group(2).splitlines() if l.strip() and not re.match(r'^\[\d+\]:$', l.strip())]
        sols = []
        for k in range(0, len(vals) - len(vs) + 1, len(vs)):
            sols.append({vs[j]: parse_complex(vals[k + j]) for j in range(len(vs))})
        res[m.group(1)]['sols'] = sols
    return res
def parse_complex(s):
    s = s.strip().replace(' ', '')
    s2 = re.sub(r'(?<![0-9a-zA-Z_.])i(?![0-9a-zA-Z_])', '1j', s)
    return complex(eval(s2, {'__builtins__': {}}, {}))
def real_up_to_torus(lay, sol):
    def am(z): return cmath.phase(z) % cmath.pi
    def same(a, b): return min(abs(a - b), cmath.pi - abs(a - b)) < 1e-15
    thu = am(sol['u1']); thv = am(sol['v1'])
    for j in range(1, lay.c + 1):
        if not same(am(sol['u%d' % j]), thu) or not same(am(sol['v%d' % j]), thv): return False
    if not same((thv - thu) % cmath.pi, 0.0): return False
    if not same(am(sol['t2']), 0.0): return False
    if lay.b_axis == 3 and not same(thu, 0.0): return False
    if lay.b_axis == 2 and not same(thv, 0.0): return False
    return True
if __name__ == '__main__':
    fn = sys.argv[1]; start, end = int(sys.argv[2]), int(sys.argv[3])
    data = json.load(open(fn))
    lay = Layout(4, data['cp'], data['b_axis'], data['bp']); vs = variables(lay)
    keys = sorted(data['structures'])
    resfile = fn.replace('.json', '_results.json')
    try: results = json.load(open(resfile))
    except Exception: results = {}
    script = 'LIB "solve.lib";\n'; tags = []
    for idx in range(start, min(end, len(keys))):
        key = keys[idx]
        if key in results: continue
        eqs, ineqs = build(lay, data['structures'][key])
        if eqs is None: results[key] = {'status': 'contradiction'}; continue
        script += singular_script(lay, eqs, 'S%d' % idx); tags.append((idx, key, ineqs))
    script += 'quit;\n'; open('/tmp/gen.sing', 'w').write(script)
    t0 = time.time()
    try:
        out = subprocess.run(['Singular', '-q', '/tmp/gen.sing'], capture_output=True, text=True, timeout=90).stdout
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or ''); print("TIMEOUT")
    parsed = parse_output(out, vs)
    for idx, key, ineqs in tags:
        tag = 'S%d' % idx
        if tag not in parsed: continue
        r = parsed[tag]; good = []
        for sol in r['sols']:
            env = dict(sol)
            if any(abs(eval(poly_str(q).replace('^', '**'), {}, env)) < 1e-12 for q in ineqs): continue
            good.append({'real': real_up_to_torus(lay, sol), 'sol': {k: [v.real, v.imag] for k, v in sol.items()}})
        results[key] = {'status': 'ok', 'dim0': r['dim0'], 'dim': r['dim'], 'vdim': r['vdim'], 'nsols': len(r['sols']),
                        'nondeg': len(good), 'nreal': sum(1 for g in good if g['real']), 'sols': good}
    json.dump(results, open(resfile, 'w'))
    done = len(results)
    print("%s: solved %d-%d in %.1fs; results %d/%d" % (fn, start, end, time.time() - t0, done, len(keys)))
    for idx, key, ineqs in tags:
        r = results.get(key)
        if r and r['status'] == 'ok' and r['dim'] >= 0:
            print("  ", idx, "dim", r['dim'], "vdim", r['vdim'], "nondeg", r['nondeg'], "real", r['nreal'])
