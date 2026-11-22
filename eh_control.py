import os
"""Certified two-probe dimension gap on one existing EH record."""
from fractions import Fraction as F
import re
import duckdb
from flint import acb,arb,ctx
ctx.prec=192
c=duckdb.connect(config={'threads':1,'memory_limit':'128MB'})
raw,ph,N,dual=c.execute('select s_matrix,twist_phases,fusion,duals from read_parquet(?) where id=?',[os.environ.get('MTC_ATLAS_MODULAR','/home/agents/data/mtc/anyon-atlas/data/modular.parquet'),'B22-00332']).fetchone()
token=re.compile(r'([+-]?)([0-9]+(?:/[0-9]+)?)?\*?(E\((-?[0-9]+/[0-9]+)\))?')
def parse(v):
 v=v.replace(' ','');ans=acb(0);pos=0
 while pos<len(v):
  m=token.match(v,pos);assert m and m.end()>pos
  sign,q,root,t=m.groups();q=F(q or 1)*(-1 if sign=='-' else 1)
  val=arb(q.numerator)/q.denominator
  if t:
   t=F(t);val*=acb(0,2*arb.pi()*t.numerator/t.denominator).exp()
  ans+=val;pos=m.end()
 return ans
H=[[parse(v) for v in row]for row in raw];d=[v.real for v in H[0]]
assert all(v.imag.contains(0) for row in H for v in row)
assert all(v>48 for v in d[1:])
assert 48<d[1]<arb(97)/2 and 49<d[5]<50
assert all(114<d[a]<115 for a in [6,7,8])
ph=[F(t)%1 for t in ph];bosons=[a for a in range(1,22) if ph[a]==0]
assert ph[1]==F(1,13) and ph[2]==F(5,13)
assert bosons==[5,6,7,8,21] and all(dual[a]==a for a in bosons)
Delta=sum(v*v for v in d)
def verlinde(i,j,k):
 v=sum(H[i][z]*H[j][z]*H[dual[k]][z]/d[z]for z in range(22))/Delta
 n=N[i][j][k]
 assert v.imag.contains(0) and v.real.contains(n)
 assert n-arb(1)/4<v.real<n+arb(1)/4
 return n
for a in bosons:
 Ds=[]
 for x in [1,2]:
  for b in range(22):verlinde(a,x,b)
  channels=[(b,N[a][x][b]) for b in range(22)if N[a][x][b] and ph[b]==ph[x]]
  assert channels==([] if a==21 else [(x,1)])
  Ds.append(sum(N[a][x][b]*d[b] for b in range(22)if ph[b]==ph[x]))
 if a!=21:
  trace=Ds[0]*Ds[1]/d[a]
  assert trace<48
  hopf=H[a][3]
  assert abs(hopf)<129 and d[a]*d[3]>2352
 assert N[a][1][2]==({5:0,6:1,7:1,8:1,21:2}[a])
print('PASS: 220 EH coefficients independently isolated; all two-probe matched channels.')
print('PASS: every nonunit simple dimension >48; all four nonzero trace bounds <48.')
print('PASS: probe 3 is a strict noncentral meridian for bosons 5,6,7,8.')
for b in range(22):verlinde(1,2,b);verlinde(3,4,b)
vac=sum(N[1][2][b]*N[3][4][dual[b]]for b in range(22))
assert vac==17
print('PASS: 44 further fusion coefficients independently isolated; vacuum(1,2,3,4)=17.')
v=[1]+[0]*21
for x in range(22):v=[sum(v[a]*N[a][x][b]for a in range(22))for b in range(22)]
assert v[0]>0
print('PASS: the full-colour word has nonzero vacuum, using the supplied fusion tensor.')
print('No symbols, large braid matrices or new modular data were generated.')
