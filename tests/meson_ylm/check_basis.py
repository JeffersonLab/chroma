"""Usage: python check_basis.py basis.txt [corr_graph.inspect.xml]. Requires numpy."""
import sys, itertools, re, xml.etree.ElementTree as ET
import numpy as np
basis={}; eta={}
for line in open(sys.argv[1]):
 a,b=line.split(':')
 if a[0]=='K': key=tuple(map(int,a[1:].split())); basis[key]={};eta[key]=int(b)
 else: basis[key][tuple(map(int,a[1:].split()))]=complex(*map(float,b.split()))
# Orthonormality of the complete Cartesian-to-coupled transformation at each order.
for n in range(4):
 paths=list(itertools.product((1,2,3),repeat=n));ks=[k for k in basis if k[0]==n]
 B=np.array([[basis[k].get(p,0) for k in ks] for p in paths])
 assert np.max(abs(B.conj().T@B-np.eye(3**n)))<1e-13
 if n:
  for k in ks:
   km=(*k[:-1],-k[-1]);factor=eta[k]*(-1)**k[-1]
   for p in paths:
    assert abs((-1)**n*basis[k].get(p[::-1],0).conjugate()-factor*basis[km].get(p,0))<1e-13
print('PASS: all 40 C++ components orthonormal; ordered adjoint identities')
# Noncommuting SU(3) links, periodic 3^3 lattice, five orthonormal vector columns.
rng=np.random.default_rng(7721);L=3;coords=list(itertools.product(range(L),repeat=3));index={x:i for i,x in enumerate(coords)};N=len(coords)*3
T=[]
for mu in range(3):
 t=np.zeros((N,N),complex)
 for x in coords:
  q,r=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)));q[:,0]/=np.linalg.det(q)
  y=list(x);y[mu]=(y[mu]+1)%L;i=index[x]*3;j=index[tuple(y)]*3;t[i:i+3,j:j+3]=q
 T.append(t)
V=np.linalg.qr(rng.normal(size=(N,5))+1j*rng.normal(size=(N,5)))[0]
def phase(p):return np.repeat(np.exp(2j*np.pi*np.array(coords)@np.array(p)/L),3)
def evaluate(p,q=(0,0,0)):
 W=phase(q)[:,None]*V;P=phase(p);A=[]
 for mu,t in enumerate(T):
  z=np.exp(2j*np.pi*p[mu]/L);a=(1+z.conjugate())*t.conj().T-(1+z)*t;A.append(a)
  D=t-t.conj().T
  explicit=(D@W).conj().T@(P[:,None]*W)-W.conj().T@(P[:,None]*(D@W))
  assert np.max(abs(explicit-W.conj().T@(P[:,None]*(a@W))))<1e-13
 paths={():W};raw={}
 for n in range(4):
  for path in itertools.product((1,2,3),repeat=n):
   if path:paths[path]=A[path[-1]-1]@paths[path[:-1]]
   raw[path]=W.conj().T@(P[:,None]*paths[path])
 return {k:sum(c*raw[p] for p,c in terms.items()) for k,terms in basis.items()}
worst=0
for q in [(0,0,0),(1,0,-1)]:
 for p in [(0,0,0),(1,0,0),(1,1,1),(1,-1,0)]:
  a=evaluate(p,q);b=evaluate(tuple(-x for x in p),q)
  for k,v in a.items():
   km=(*k[:-1],-k[-1]) if k[0] else k
   sign=eta[k]*(-1)**k[-1] if k[0] else 1
   worst=max(worst,float(np.max(abs(v.conj().T-sign*b[km]))))
assert worst<1e-11
print('PASS: rectangular-column SU(3) test, left-to-right conversion and opposite-momentum reconstruction; max error',worst)
if len(sys.argv)>2:
 records=[];pairs={};maxres=0;terms=0
 for e in ET.parse(sys.argv[2]).getroot().find('VertexCoeffMap'):
  k=e.find('Key');name=k.findtext('name');n,j13,J=re.search(r'xD(\d+)(?:_J13(\d+))?_J(\d+)__',name).groups();n=int(n);J=int(J);j13=int(j13) if j13 else None
  paths=list(itertools.product((1,2,3),repeat=n));keys=[x for x in basis if x[0]==n and (n==1 or x[-2]==J) and (n!=3 or x[1]==j13)]
  B=np.array([[basis[x].get(p,0) for x in keys] for p in paths]);d={}
  for f in e.findall('Val/elem'):
   for t in f.findall('Val/elem'):
    ds=[tuple(map(int,(x.text or '').split())) for x in t.findall('Key/derivs/elem/deriv')];assert len(ds)==2 and not ds[0] and len(ds[1])==n
    vals={tuple(map(int,s.findtext('Key/diracs').split())):complex(float(s.findtext('Val/re')),float(s.findtext('Val/im'))) for s in t.findall('Val/elem')}
    assert set(vals)=={(i,i) for i in range(4)} and len(set(vals.values()))==1
    d[ds[1]]=vals[0,0];terms+=len(vals)
  c=np.array([d.get(p,0) for p in paths]);w=B.conj().T@c;maxres=max(maxres,float(np.max(abs(c-B@w))));assert abs(np.linalg.norm(w)-1)<1e-12
  identity=tuple((x.tag,ET.tostring(x,encoding='unicode').strip()) for x in k if x.tag!='creation_op');pairs.setdefault(identity,{})[k.findtext('creation_op')]=d
 for pair in pairs.values():
  assert set(pair)=={'true','false'}
  expect={p[::-1]:(-1)**len(p)*v.conjugate() for p,v in pair['true'].items()}
  assert set(expect)==set(pair['false'])
  assert max(abs(expect[p]-pair['false'][p]) for p in expect)<1e-13
 assert maxres<1e-12
 print('PASS:',2*len(pairs),'Redstar vertices,',terms,'spin coefficients; C++ basis residual',maxres)
