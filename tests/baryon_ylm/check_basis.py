"""Independent baryon CG/circular and finite-momentum SU(3) checks; numpy required.
Usage: python check_basis.py basis.txt
"""
import itertools,sys
import numpy as np
basis={}
for line in open(sys.argv[1]):
 if line.startswith('K'):
  key=tuple(map(int,line[1:].split()));basis[key]={}
 else:
  p,w=line.split(':');basis[key][tuple(tuple(map(int,q.split())) for q in p.split('|'))]=complex(*map(float,w.split()))
assert len(basis)==22
# Independent 1 x 1 CG table, physical integer labels (Condon--Shortley).
def cg(L,M,a,b):
 if a+b!=M:return 0.
 if L==0:return (-1.)**(1-a)/np.sqrt(3)
 if L==1:return (1 if a>b else -1 if a<b else 0)/np.sqrt(2)
 if abs(M)==2:return 1.
 if abs(M)==1:return 1/np.sqrt(2)
 return (2 if a==b==0 else 1)/np.sqrt(6)
circ={1:{1:.5j,2:-.5},0:{3:-1j/np.sqrt(2)},-1:{1:-.5j,2:-.5}}
for key,terms in basis.items():
 expected={}
 if key==(0,):expected[((),(),())]=1
 elif key[0]==3:expected={((),(),(d,)):w for d,w in circ[key[1]].items()}
 else:
  f,L,M=key
  for a,b in itertools.product(range(-1,2),repeat=2):
   for x,wx in circ[a].items():
    for y,wy in circ[b].items():
     p=((),(),(y,x)) if f==33 else ((),(x,),(y,))
     expected[p]=expected.get(p,0)+cg(L,M,a,b)*wx*wy
 assert max(abs(expected.get(p,0)-terms.get(p,0)) for p in set(expected)|set(terms))<1e-13
for f,n in [(0,0),(3,1),(33,2),(23,2)]:
 ks=[k for k in basis if k[0]==f];ps=sorted({p for k in ks for p in basis[k]})
 A=np.array([[basis[k].get(p,0) for k in ks] for p in ps])
 assert np.max(abs(A.conj().T@A-np.eye(len(ks))*2.**(-n)))<1e-13
rng=np.random.default_rng(927);L=3;coords=list(itertools.product(range(L),repeat=3));idx={x:i for i,x in enumerate(coords)};N=3*len(coords);nv=4
T=[]
for mu in range(3):
 t=np.zeros((N,N),complex)
 for x in coords:
  u,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)));u[:,0]/=np.linalg.det(u)
  y=list(x);y[mu]=(y[mu]+1)%L;i=idx[x]*3;j=idx[tuple(y)]*3;t[i:i+3,j:j+3]=u
 T.append(t)
D=[t-t.conj().T for t in T];d={m:sum(w*D[x-1] for x,w in ws.items()) for m,ws in circ.items()}
V=np.linalg.qr(rng.normal(size=(N,nv))+1j*rng.normal(size=(N,nv)))[0]
eps=np.zeros((3,3,3))
for a,b,c in itertools.permutations(range(3)):eps[a,b,c]=((b-a)*(c-a)*(c-b))//2
worst=0.;permworst=0.
for q in [(0,0,0),(1,0,-1)]:
 W=np.repeat(np.exp(2j*np.pi*np.array(coords)@q/L),3)[:,None]*V
 paths={():W}
 for n in [1,2]:
  for p in itertools.product(range(1,4),repeat=n):paths[p]=D[p[-1]-1]@paths[p[:-1]]
 for mom in [(0,0,0),(1,0,0),(-1,0,0),(1,1,1),(-1,-1,-1)]:
  phase=np.exp(2j*np.pi*np.array(coords)@mom/L)
  def C(a,b,c):return np.einsum('x,abc,xai,xbj,xck->ijk',phase,eps,a.reshape(-1,3,nv),b.reshape(-1,3,nv),c.reshape(-1,3,nv),optimize=True)
  for k,terms in basis.items():
   value=sum(w*C(*(paths[p] for p in ps)) for ps,w in terms.items())
   if k==(0,):expected=C(W,W,W)
   elif k[0]==3:expected=C(W,W,d[k[1]]@W)
   else:
    f,J,M=k;expected=np.zeros_like(value)
    for a,b in itertools.product(range(-1,2),repeat=2):
     expected+=cg(J,M,a,b)*(C(W,W,d[a]@d[b]@W) if f==33 else C(W,d[a]@W,d[b]@W))
   worst=max(worst,float(np.max(abs(value-expected))))
   partner=value.swapaxes(1,2) if k[0]==23 else value.swapaxes(0,1)
   sign=(-1)**(k[1]+1) if k[0]==23 else -1
   permworst=max(permworst,float(np.max(abs(partner-sign*value))))
assert worst<1e-12 and permworst<1e-12
print('PASS: 22 components; independent CG/circular coefficients and normalized Gram matrices')
print('PASS: finite-momentum SU(3) direct construction, max error',worst)
print('PASS: same-line and split-line permutation identities, max error',permworst)
