import os
"""Exact twist-support inputs for the proposed gauging closure.

The quantitative boson bounds are those of the fully checked T494/T496.
Only two supplied records are read; no category or symbol data generated.
"""
from fractions import Fraction as F
import duckdb
c=duckdb.connect(config={'threads':1,'memory_limit':'128MB'})
for ident,probe,q,bosons in [('B22-00332',9,65,[5,6,7,8,21]),('B16-00499',15,26,[5,10])]:
 phases,N,dual=c.execute('select twist_phases,fusion,duals from read_parquet(?) where id=?',[os.environ.get('MTC_ATLAS_MODULAR','/home/agents/data/mtc/anyon-atlas/data/modular.parquet'),ident]).fetchone()
 phases=[F(t)%1 for t in phases];support=set(phases)
 orders={t.denominator for t in phases}
 assert orders==({1,5,13}if q==65 else{1,2,13})
 assert [a for a in range(1,len(phases))if phases[a]==0]==bosons
 assert phases[1].denominator==13 and phases[probe].denominator==(5 if q==65 else 2)
 for a in range(1,len(phases)):
  if phases[a]==0:continue
  x=probe if phases[a].denominator==13 else 1
  product=(phases[a]+phases[x])%1
  assert product.denominator==q and product not in support
 print('PASS:',ident,'all nonbosons have an empty M_1 packet with probe1 or',probe)
 if q==26:
  assert [(b,N[dual[5]][1][b])for b in range(16)if N[dual[5]][1][b] and phases[b]==phases[1]]==[(1,1)]
  assert not any(N[dual[10]][1][b]for b in range(16)if phases[b]==phases[1])
  assert 19*F(8,18)**4<1
  print('PASS: F4 boson channels and exact four-copy scalar majorant.')
print('No new modular data, symbols or braid matrices produced.')
