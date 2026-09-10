import json, itertools, numpy as np, sys
from scipy.optimize import least_squares
np.random.seed(1)
types={}
for fn in ['types_22_11_13_17.json','types_22_19.json']:
    for k,v in json.load(open(fn)).items(): types.setdefault(k,[]).extend(v)
key=[k for k,v in types.items() if len({o['p'] for o in v})>1][0]
occ=[o for o in types[key] if o['p']==13][0]
p=13; B=occ['form']; P=[tuple(x) for x in occ['points']]
print("lifting type found at primes",sorted({o['p'] for o in types[key]}),"from p=13, form",B)
n=len(P); assert n==22
import re
if B.startswith('diag'):
    d=[int(x) for x in re.findall(r'-?\d+',B)]; G=[[d[0],0,0],[0,d[1],0],[0,0,d[2]]]
else:
    c=int(re.findall(r'\d+',B)[1])//2; G=[[1,0,0],[0,0,c],[0,c,0]]
def bil(u,v): return sum(G[i][j]*u[i]*v[j] for i in range(3) for j in range(3))%p
inc=[[bil(P[i],P[j])==0 for j in range(n)] for i in range(n)]
def kind(q):
    z=[c==0 for c in q]
    if z==[True,False,False]: return 'l1'
    if z==[False,True,False]: return 'l2'
    if z==[False,False,True]: return 'l3'
    if sum(z)==0: return 'reg'
    return 'coord'
kinds=[kind(q) for q in P]; print("kinds:",{k:kinds.count(k) for k in set(kinds)})
assert 'coord' not in kinds
V=[(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
def norm(q):
    for c in q:
        if c:
            ic=pow(c,-1,p); return tuple((ic*x)%p for x in q)
orb_of={}; orbits=[]
for i,q in enumerate(P):
    if i in orb_of: continue
    o=[]
    for v in V:
        r=norm(tuple(s*c for s,c in zip(v,q)))
        j=P.index(r);
        if j not in o: o.append(j)
    for j in o: orb_of[j]=len(orbits)
    orbits.append(o)
print("orbits:",[len(o) for o in orbits])
params=[]; templ={}
for oi,o in enumerate(orbits):
    k=kinds[o[0]]
    if k=='reg': idx=(len(params),len(params)+1); params+=['u%d'%oi,'v%d'%oi]
    else: idx=(len(params),); params+=['t%d'%oi]
    templ[oi]=(k,idx)
signs={}
for oi,o in enumerate(orbits):
    rep=P[o[0]]
    for j in o:
        for v in V:
            if norm(tuple(s*c for s,c in zip(v,rep)))==P[j]: signs[j]=v; break
def coords(x, j):
    oi=orb_of[j]; k,idx=templ[oi]; v=signs[j]
    if k=='l1': q=np.array([0,1,x[idx[0]]],dtype=complex)
    elif k=='l2': q=np.array([1,0,x[idx[0]]],dtype=complex)
    elif k=='l3': q=np.array([1,x[idx[0]],0],dtype=complex)
    else: q=np.array([1,x[idx[0]],x[idx[1]]],dtype=complex)
    q=q*np.array(v)
    for c in q:
        if abs(c)>0: q=q/c; break
    return q
pairs=[(i,j) for i in range(n) for j in range(i,n) if inc[i][j]]
nonpairs=[(i,j) for i in range(n) for j in range(i,n) if not inc[i][j]]
def resid(xr):
    x=xr[:len(params)]+1j*xr[len(params):]
    Q=[coords(x,j) for j in range(n)]
    r=[Q[i]@Q[j] for i,j in pairs]
    r=np.array(r); return np.concatenate([r.real,r.imag])
sols=[]
for trial in range(int(sys.argv[1]) if len(sys.argv)>1 else 60):
    x0=np.random.randn(2*len(params))*1.5
    res=least_squares(resid,x0,xtol=1e-14,ftol=1e-14,gtol=1e-14,max_nfev=1500)
    if res.cost<1e-20:
        x=res.x[:len(params)]+1j*res.x[len(params):]
        Q=[coords(x,j) for j in range(n)]
        if min(abs(Q[i]@Q[j]) for i,j in nonpairs)<1e-6: continue
        if min(np.linalg.norm(np.cross(Q[i],Q[j])) for i,j in itertools.combinations(range(n),2))<1e-6: continue
        if not any(np.allclose(x,s,atol=1e-7) for s in sols): sols.append(x)
print("nondegenerate solutions found:",len(sols))
for s in sols[:12]:
    print("  ",[ "%s=%.6f%+.6fi"%(pn,z.real,z.imag) for pn,z in zip(params,s)])
json.dump({'params':params,'sols':[[ [z.real,z.imag] for z in s] for s in sols],'templ':{str(k):v for k,v in templ.items()},
           'orbits':orbits,'signs':{str(k):v for k,v in signs.items()},'kinds':kinds,'inc':inc},open('lift22.json','w'))
