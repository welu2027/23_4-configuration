import json, collections, itertools
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher

def levi(m):
    G = nx.Graph(); k = len(m)
    for i in range(k): G.add_node(("p", i), side=0); G.add_node(("l", i), side=1)
    for i in range(k):
        for j in range(k):
            if m[i][j]: G.add_edge(("p", i), ("l", j))
    return G

def report(name, m):
    k = len(m)
    assert all(sum(r) == 4 for r in m), name + ": some point is not on 4 lines"
    assert all(sum(m[i][j] for i in range(k)) == 4 for j in range(k)), name + ": some line lacks 4 points"
    assert all(sum(a and b for a, b in zip(m[i], m[j])) <= 1
               for i, j in itertools.combinations(range(k), 2)), name + ": two points on two lines"
    G = levi(m); sides = lambda a, b: a["side"] == b["side"]
    apart = sum(1 for _ in GraphMatcher(G, G, node_match=sides).isomorphisms_iter())
    alld  = sum(1 for _ in GraphMatcher(G, G).isomorphisms_iter())
    print("%s: geometric (%d_4), Levi girth %d" % (name, k, nx.girth(G)))
    print("   Aut(Levi) with points and lines kept apart: %d" % apart)
    print("   Aut(Levi) with dualities allowed:           %d" % alld)
    print("   self-dual: %s   (hence self-polar: %s)" % (alld > apart, "possible" if alld > apart else "no"))
    return apart, alld

from fractions import Fraction as F
from itertools import combinations
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher

D=17
def mk(a,b=0): return (F(a),F(b))
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def sub(x,y): return (x[0]-y[0],x[1]-y[1])
def mul(x,y): return (x[0]*y[0]+D*x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def iszero(x): return x[0]==0 and x[1]==0
def inv(x):
    n=x[0]*x[0]-D*x[1]*x[1]; return (x[0]/n,-x[1]/n)
Z=mk(0); ONE=mk(1)
w=(F(-7,2),F(3,2))
assert iszero(add(add(mul(w,w),mul(mk(7),w)),mk(-26)))
def lin(a,b): return add(mul(mk(a),w),mk(b))

L=[(mk(1),mk(0),mk(0)),(mk(0),mk(1),mk(0)),(mk(0),mk(0),mk(1)),(mk(1),mk(1),mk(1)),
   (mk(24),lin(-5,-13),mk(0)),(mk(24),lin(5,13),lin(24,0)),(mk(1),mk(0),w),(mk(2),mk(0),w),
   (mk(24),lin(-5,-13),lin(-4,52)),(mk(24),lin(5,13),lin(28,-52)),
   (mk(6),lin(-1,13),lin(-1,13)),(mk(24),lin(-5,-13),lin(16,104)),
   (mk(48),lin(1,65),lin(24,0)),(mk(24),lin(5,13),lin(-32,104)),
   (mk(18),lin(-1,13),lin(4,26)),(mk(12),lin(-1,13),mk(0)),
   (mk(96),lin(1,65),lin(56,-104)),(mk(48),lin(1,65),lin(-8,104)),
   (mk(48),lin(1,65),lin(20,52)),(mk(39),lin(-1,52),lin(-1,52)),
   (mk(4),lin(1,13),lin(4,0)),(mk(24),lin(1,26),lin(12,0))]
assert len(L)==22
def cross(u,v):
    return (sub(mul(u[1],v[2]),mul(u[2],v[1])),sub(mul(u[2],v[0]),mul(u[0],v[2])),sub(mul(u[0],v[1]),mul(u[1],v[0])))
def norm(p):
    for c in p:
        if not iszero(c):
            ic=inv(c); return tuple(mul(ic,x) for x in p)
def dot(p,l): return add(add(mul(p[0],l[0]),mul(p[1],l[1])),mul(p[2],l[2]))
assert all(any(not iszero(t) for t in cross(L[i],L[j])) for i,j in combinations(range(22),2))
pts={}
for i,j in combinations(range(22),2):
    p=norm(cross(L[i],L[j])); pts.setdefault(p,set()).update({i,j})
mult={}
for p,s in pts.items(): mult.setdefault(len(s),0); mult[len(s)]+=1
print("intersection points by multiplicity:",dict(sorted(mult.items())))
Q=[p for p,s in pts.items() if len(s)==4]
print("quadruple points:",len(Q))
inc=[[iszero(dot(p,l)) for l in L] for p in Q]
print("each line has 4 quadruple points:", all(sum(inc[i][j] for i in range(len(Q)))==4 for j in range(22)))
print("each quadruple point on 4 lines:", all(sum(r)==4 for r in inc))
print("no two points on two lines:", all(sum(a and b for a,b in zip(inc[i],inc[j]))<=1 for i,j in combinations(range(len(Q)),2)))

INC22 = inc
report("Cuntz (22_4)", INC22)
json.dump([[int(x) for x in r] for r in INC22], open('cuntz22_inc.json','w'))

from fractions import Fraction as F
from itertools import combinations
MOD=[F(7),F(1),F(-3)]
def C(*c): return [F(x) for x in c]+[F(0)]*(3-len(c))
Z=C(0); ONE=C(1); zz=C(0,1)
def add(x,y): return [a+b for a,b in zip(x,y)]
def sub(x,y): return [a-b for a,b in zip(x,y)]
def mul(x,y):
    r=[F(0)]*5
    for i in range(3):
        for j in range(3): r[i+j]+=x[i]*y[j]
    for k in (4,3):
        if r[k]:
            c=r[k]; r[k]=F(0)
            for i in range(3): r[k-3+i]+=c*MOD[i]
    return r[:3]
def iszero(x): return all(a==0 for a in x)
def inv(x):
    M=[[F(0)]*3 for _ in range(3)]
    for j in range(3):
        e=[F(0)]*3; e[j]=F(1); p=mul(x,e)
        for i in range(3): M[i][j]=p[i]
    A=[row[:]+[F(1) if i==0 else F(0)] for i,row in enumerate(M)]
    for c in range(3):
        piv=next(r for r in range(c,3) if A[r][c]!=0); A[c],A[piv]=A[piv],A[c]
        pv=A[c][c]; A[c]=[v/pv for v in A[c]]
        for r in range(3):
            if r!=c and A[r][c]!=0:
                f=A[r][c]; A[r]=[a-f*b for a,b in zip(A[r],A[c])]
    return [A[i][3] for i in range(3)]
def P(*coeffs):
    return C(*coeffs)

L=[
 (C(1),C(0),C(0)), (C(0),C(1),C(0)), (C(0),C(0),C(1)), (C(1),C(1),C(1)),
 (C(1), P(0,-2,-1), P(0,1)), (P(0,1), P(0,-2,-1), P(0,1)), (C(-1), P(0,-1), P(0,-1)),
 (P(0,2,2), P(0,-14,-6), P(21,-4,-5)),
 (P(-14,12,10), P(-98,-56,-6), P(-28,24,16)),
 (P(-14,0,6), P(14,-12,-2), P(56,-20,-16)),
 (P(-196,56,68), P(-84,72,20), P(-280,100,84)),
 (P(-28,-60,-24), C(0), P(14,-40,-26)),
 (P(784,-112,-256), P(336,-624,-352), P(392,-224,-264)),
 (C(0), P(28,4,-16), P(-14,12,-2)),
 (P(-196,56,68), P(-84,72,20), P(-84,72,20)),
 (P(3696,-256,-1136), P(7952,-2560,-3152), P(4928,-528,-1840)),
 (P(1792,-1760,-608), P(-8064,-1824,1120), P(6048,-4064,-2624)),
 (C(0), P(8064,1824,-1120), P(7056,-1120,-1872)),
 (P(48832,-14976,-12864), P(90048,-35968,-23616), P(88256,-19200,-27584)),
 (P(61376,37888,4288), P(-157248,-170752,-44736), P(48608,19712,-2656)),
 (P(-224,136,8), P(1176,-392,-304), P(1652,-632,-412)),
 (P(2800,-608,-784), P(-1232,-288,-272), P(3696,-256,-1136)),
 (P(-75264,30464,11264), P(637952,-190208,-193536), P(307328,-57344,-123776)),
 (P(231616,-13056,-65984), P(37184,-91904,-55360), P(-20160,55808,448)),
 (P(-155904,-31232,8192), P(-57344,147712,55808), P(17024,-36096,-54912)),
 (P(-1732864,367104,627968), P(-9105152,3188224,2495232), P(-5465600,1928704,1453568)),
]
def cross(u,v):
    return (sub(mul(u[1],v[2]),mul(u[2],v[1])), sub(mul(u[2],v[0]),mul(u[0],v[2])),
            sub(mul(u[0],v[1]),mul(u[1],v[0])))
def norm(p):
    for c in p:
        if not iszero(c):
            ic=inv(c); return tuple(tuple(mul(ic,x)) for x in p)
def dot(p,l): return add(add(mul(p[0],l[0]),mul(p[1],l[1])),mul(p[2],l[2]))

pts = collections.defaultdict(set)
for i, j in itertools.combinations(range(26), 2):
    pts[norm(cross(L[i], L[j]))].update({i, j})
Q = [p for p, v in pts.items() if len(v) == 4]
assert len(Q) == 26
INC26 = [[iszero(dot(p, l)) for l in L] for p in Q]
report("Cuntz (26_4)", INC26)
json.dump([[int(x) for x in r] for r in INC26], open('cuntz26_inc.json','w'))
