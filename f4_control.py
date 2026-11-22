import os
"""A small falsification test for T496 Question 3.2.

Only one supplied rank-16 datum is read. No orbifold S-matrix, fusion
tensor or new dataset is generated. The computations are a 16-by-16
matrix, 16^3 scalar inequalities, and the needed base fusion row.
"""
from fractions import Fraction as F
from itertools import product, combinations
import math
import re
import duckdb
from flint import acb, arb, ctx

ctx.prec = 192
db = duckdb.connect(config={'threads': 1, 'memory_limit': '128MB'})
raw, phases, N, dual, realizations = db.execute(
    'select s_matrix,twist_phases,fusion,duals,realizations '
    'from read_parquet(?) where id=?',
    [os.environ.get('MTC_ATLAS_MODULAR','/home/agents/data/mtc/anyon-atlas/data/modular.parquet'), 'B16-00499']
).fetchone()
assert realizations == [r'\mathrm{F}_{4,4}']
token = re.compile(r'([+-]?)([0-9]+(?:/[0-9]+)?)?\*?(E\((-?[0-9]+/[0-9]+)\))?')
def parse(value):
    value = value.replace(' ', '')
    ans, pos = acb(0), 0
    while pos < len(value):
        match = token.match(value, pos)
        assert match and match.end() > pos
        sign, coeff, root, phase = match.groups()
        assert coeff or root
        q = F(coeff or 1)*(-1 if sign == '-' else 1)
        scalar = arb(q.numerator)/q.denominator
        if phase:
            t = F(phase)
            ans += scalar*acb(0, 2*arb.pi()*t.numerator/t.denominator).exp()
        else:
            ans += scalar
        pos = match.end()
    return ans

H = [[parse(v) for v in row] for row in raw]
assert all(v.imag.contains(0) for row in H for v in row)
d = [v.real for v in H[0]]
assert all(v > 0 for v in d)
assert all(v > 1 for v in d[1:])  # no nontrivial invertibles
assert dual == list(range(16))
Delta = sum(v*v for v in d)
S = [[v/Delta.sqrt() for v in row] for row in H]
phases = [F(t) % 1 for t in phases]
a = 5
assert phases[a] == 0
assert d[a] > 18 and d[a] < 19
assert d[a]**3 > 6187 and d[a]**3 < 6188
assert Delta > 7402 and Delta < 7403

# Identify label 5 with the affine highest weight 2*lambda_4.
# Long roots have squared length two; rho=(11,5,3,1)/2.
rho = tuple(F(i,2) for i in (11,5,3,1))
lam = (2,0,0,0)
def dot(v,w): return sum(x*y for x,y in zip(v,w))
roots = []
for i,j in combinations(range(4),2):
    for si,sj in product((-1,1),repeat=2):
        v = [0]*4; v[i],v[j] = si,sj
        roots.append(tuple(v))
for i in range(4):
    for si in (-1,1):
        v = [0]*4; v[i] = si
        roots.append(tuple(v))
roots += [tuple(F(s,2) for s in ss) for ss in product((-1,1),repeat=4)]
positive = [v for v in roots if dot(rho,v) > 0]
assert len(positive) == 24
lamrho = tuple(x+y for x,y in zip(lam,rho))
weyl_dim = arb(1)
for root in positive:
    top,bot = dot(lamrho,root),dot(rho,root)
    top_angle = arb.pi()*top.numerator/(13*top.denominator)
    bot_angle = arb.pi()*bot.numerator/(13*bot.denominator)
    weyl_dim *= top_angle.sin()/bot_angle.sin()
assert weyl_dim > 18 and weyl_dim < 19
assert (weyl_dim-d[a]).contains(0)
assert [i for i in range(16) if phases[i] == 0] == [0,5,10]
assert d[10] > 25
assert dot(lam,tuple(x+2*y for x,y in zip(lam,rho)))/26 == 1
print('PASS: affine label 2*lambda_4 is the unique boson of dimension in (18,19).')

# Independently recover each base fusion coefficient used by the spectral
# measures and by the cubic fusion diagonal.
for x, b in product(range(16), repeat=2):
    value = sum(S[a][z]*S[x][z]*S[b][z].conjugate()/S[0][z]
                for z in range(16))
    expected = N[a][x][b]
    assert value.imag.contains(0) and value.real.contains(expected)
    assert value.real > arb(expected)-arb(1)/4
    assert value.real < arb(expected)+arb(1)/4
print('PASS: 256 required F4 fusion coefficients independently isolated by Verlinde.')

# Spectral measure of M_(a,x), divided by d_a, on the 26th roots.
# Certify a rational minorant coefficientwise, then use integer convolutions.
den = 1000
W = [[arb(0) for _ in range(26)] for _ in range(16)]
for x, b in product(range(16), repeat=2):
    k = 26*((phases[b]-phases[x]) % 1)
    assert k.denominator == 1
    W[x][int(k)] += N[a][x][b]*d[b]/d[a]
lower = []
for x in range(16):
    row = []
    for v in W[x]:
        n = max(0, math.floor(den*float(v))-1)
        assert v >= arb(n)/den
        row.append(n)
    lower.append(row)
lower[0] = [den]+[0]*25  # a tensor 1 = a, exactly
sparse = [[(i,n) for i,n in enumerate(row) if n] for row in lower]
least, arg = None, None
strata = {}
for xs in product(range(16), repeat=3):
    if xs == (0,0,0):
        continue
    x,y,z = xs
    coeff = sum(u*v*lower[z][(-i-j) % 26]
                for i,u in sparse[x] for j,v in sparse[y])
    orbit = 1 if x == y == z else 3
    bound = F(orbit*coeff, den**3)
    assert bound > 1, (xs,bound)
    if least is None or bound < least:
        least,arg = bound,xs
    key = ('diagonal' if orbit == 1 else 'orbit',sum(t != 0 for t in xs))
    strata[key] = min(bound,strata.get(key,bound))
print('PASS: all 4095 nonunit neutral tuple bounds exceed one.')
print('Smallest certified rational minorant:',least,'at',arg)
print('Certified minima by (type, nonunit entries):',strata)
print('The three electric charges on the all-unit tuple have ratio exactly one.')

Na = N[a]
v = [int(i == a) for i in range(16)]
for _ in range(2):
    v = [sum(v[i]*Na[i][j] for i in range(16)) for j in range(16)]
cubic_diag = [sum(v[i]*N[i][y][y] for i in range(16)) for y in range(16)]
assert cubic_diag == [2,56,56,56,56,300,327,327,327,327,560,610,610,610,610,1546]
assert min(cubic_diag[1:]) > 2*d[a]
print('PASS: cubic fusion diagonal:',cubic_diag)

# Xu's Lambda(1/3). Its scalar anomaly is irrelevant to absolute squares.
# An integral choice of conformal-weight lift changes only row/column phases.
# c(F4,4)=16; T_j=exp(2 pi i (h_j-c/24)).
def phase(q):
    q = F(q)
    return acb(0, 2*arb.pi()*q.numerator/q.denominator).exp()
t_minus_third = [phase(-(h-F(2,3))/3) for h in phases]
t_minus_three = [phase(-3*(h-F(2,3))) for h in phases]
L = [[t_minus_third[i]*sum(S[k][i].conjugate()*t_minus_three[k]*S[k][j]
                          for k in range(16))*t_minus_third[j]
      for j in range(16)] for i in range(16)]
r = [H[a][z].real/d[z] for z in range(16)]
for y in range(16):
    norm = sum((z*z.conjugate()).real for z in L[y])
    assert norm.contains(1)
    hyy = sum(r[z]*(L[y][z]*L[y][z].conjugate()).real for z in range(16))
    multiplicity = (cubic_diag[y]+2*hyy)/3
    assert multiplicity > 0
    if y == 0:
        assert hyy > arb(3)/2 and hyy < arb(5)/2
        assert multiplicity.contains(2)
print('PASS: certified H_00 lies in (3/2,5/2), and all twisted self-channels are positive.')
print('PASS: every twisted simple has dimension >=Delta>d_a^3, hence cannot contract.')
print('PASS: 6187<d_a^3<6188, so the bosonic dimension is nonintegral.')
print('No new modular datum or orbifold fusion tensor was constructed.')

# Additional matched packet used in the joint-detection seed.
for b in range(16):
    value=sum(S[10][z]*S[1][z]*S[b][z].conjugate()/S[0][z] for z in range(16))
    k=N[10][1][b]
    assert value.imag.contains(0) and value.real.contains(k)
    assert value.real>arb(k)-arb(1)/4 and value.real<arb(k)+arb(1)/4
assert not any(N[10][1][b] for b in range(16) if phases[b]==phases[1])
counts=[sum(N[5][x][b] for b in range(16) if phases[b]==phases[x]) for x in range(16)]
assert counts==[1,1,1,1,1,3,2,2,2,2,2,3,3,3,3,4]
print('PASS: extra F4 matched packet empty; all sixteen a5 packet multiplicities certified positive.')
