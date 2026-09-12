from fractions import Fraction as F
from itertools import combinations, product
D=3
def mk(a,b=0): return (F(a),F(b))
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def sub(x,y): return (x[0]-y[0],x[1]-y[1])
def mul(x,y): return (x[0]*y[0]+D*x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def inv(x):
    n=x[0]*x[0]-D*x[1]*x[1]; return (x[0]/n,-x[1]/n)
def iszero(x): return x[0]==0 and x[1]==0
s3=mk(0,1); Z=mk(0); ONE=mk(1)
def r(a,b): return mk(F(a,b))
th=lambda n,d: (F(0),F(n,d))

Preps=[ (Z, ONE, add(mk(2),s3)),
        (ONE, Z, sub(ONE, th(1,3))),
        (ONE, add(mk(-1),th(1,3)), Z),
        (ONE, ONE, ONE),
        (ONE, sub(mk(2),s3), sub(mk(2),s3)),
        (ONE, add(mk(-1),th(2,3)), ONE),
        (ONE, sub(mk(2),s3), (F(0),F(-1,3))) ]
Lreps=[ (Z, ONE, mk(-1)),
        (ONE, Z, mk(-1)),
        (ONE, sub(mk(-2),s3), Z),
        (ONE, add(r(9,2),th(5,2)), sub(r(-3,2),th(1,2))),
        (ONE, add(r(3,2),th(1,2)), add(r(-3,2),th(1,2))),
        (ONE, add(r(1,2),th(1,2)), sub(r(-3,2),th(1,2))),
        (ONE, add(r(3,2),th(1,2)), sub(r(-1,2),th(1,2))) ]

V=[(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
def act(v,p): return tuple(mul(mk(s),c) for s,c in zip(v,p))
def cross(u,v):
    return (sub(mul(u[1],v[2]),mul(u[2],v[1])),sub(mul(u[2],v[0]),mul(u[0],v[2])),
            sub(mul(u[0],v[1]),mul(u[1],v[0])))
def same(p,q): return all(iszero(t) for t in cross(p,q))
def orbit(p):
    o=[]
    for v in V:
        q=act(v,p)
        if not any(same(q,x) for x in o): o.append(q)
    return o
P=[]; sizesP=[]
for p in Preps:
    o=orbit(p); sizesP.append(len(o)); P+=o
L=[]; sizesL=[]
for l in Lreps:
    o=orbit(l); sizesL.append(len(o)); L+=o
print("point orbit sizes:",sizesP,"total",len(P))
print("line  orbit sizes:",sizesL,"total",len(L))
def dot(p,l): return add(add(mul(p[0],l[0]),mul(p[1],l[1])),mul(p[2],l[2]))
distinctP=all(not same(P[i],P[j]) for i,j in combinations(range(len(P)),2))
distinctL=all(not same(L[i],L[j]) for i,j in combinations(range(len(L)),2))
print("22 distinct points:",distinctP,"  22 distinct lines:",distinctL)
inc=[[iszero(dot(p,l)) for l in L] for p in P]
print("every point on exactly 4 lines:",all(sum(r_)==4 for r_ in inc))
print("every line has exactly 4 points:",all(sum(inc[i][j] for i in range(len(P)))==4 for j in range(len(L))))
print("no two points on two common lines:",
      all(sum(a and b for a,b in zip(inc[i],inc[j]))<=1 for i,j in combinations(range(len(P)),2)))

import networkx as nx, json, os
from networkx.algorithms.isomorphism import GraphMatcher
G=nx.Graph()
for i in range(22): G.add_node(('p',i),side=0); G.add_node(('l',i),side=1)
for i in range(22):
    for j in range(22):
        if inc[i][j]: G.add_edge(('p',i),('l',j))
sides=lambda a,b:a['side']==b['side']
print("Aut(Levi T3) kept apart:", sum(1 for _ in GraphMatcher(G,G,node_match=sides).isomorphisms_iter()),
      " with dualities:", sum(1 for _ in GraphMatcher(G,G).isomorphisms_iter()))
print("Levi girth:", nx.girth(G))
def cr(p1,p2,p3,p4):
    det=lambda a,b: sub(mul(a[0],b[1]),mul(a[1],b[0]))
    for (i,j) in ((0,1),(0,2),(1,2)):
        q=[(p[i],p[j]) for p in (p1,p2,p3,p4)]
        if all(not iszero(det(q[a],q[b])) for a,b in combinations(range(4),2)):
            return mul(mul(det(q[0],q[2]),det(q[1],q[3])), inv(mul(det(q[1],q[2]),det(q[0],q[3]))))
irr=sum(1 for j in range(22)
        for v in [cr(*[P[i] for i in range(22) if inc[i][j]])]
        if v is None or v[1]!=0)
print("collinear quadruples with irrational cross-ratio: %d of 22" % irr)
for fn in ('cuntz22_inc.json','../certificates/cuntz22_inc.json'):
    if os.path.exists(fn):
        c=json.load(open(fn)); H=nx.Graph()
        for i in range(22): H.add_node(('p',i)); H.add_node(('l',i))
        for i in range(22):
            for j in range(22):
                if c[i][j]: H.add_edge(('p',i),('l',j))
        print("T3 isomorphic to Cuntz's (22_4):", nx.is_isomorphic(G,H))
        break
else:
    print("(cuntz22_inc.json not found; run check_cuntz.py for the comparison)")
