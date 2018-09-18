#!/usr/bin/env python

from simp_sec_quant import *

# define orbitals
a = 'a' # vir_orb()
b = 'b' # vir_orb()
c = 'c' # vir_orb()
d = 'd' # vir_orb()
e = 'e' # vir_orb()
i = 'i' # occ_orb()
j = 'j' # occ_orb()
k = 'k' # occ_orb()
l = 'l' # occ_orb()
m = 'm' # occ_orb()

# define states
bHF = bra(rank=0)
bra1 = bra(rank=1, vir=[a], occ=[i])
bra2 = bra(rank=2, vir=[a,b], occ=[i,j])
kHF = ket(rank=0)
ket1 = ket(rank=1, vir=[a], occ=[i])
ket2 = ket(rank=2, vir=[a,b], occ=[i,j])

# define operators
fock = operator(rank=1, name="fock", symb="F")
phi = operator(rank=2, name="fluctuation potential", symb="\Phi")


# CIS(D) excitation energy:
lcc0 = tensor([a,i], 'L^\\text{CCS}')
rcc0 = tensor([d,l], 'R^\\text{CCS}')
tcc1 = tensor([b,j,c,k], 't^{(1)}')
c1 = commutator(phi, el_oper(b, j) * el_oper(c, k))
comm = commutator(c1, el_oper(d, l))

term1 = term(elmts=[0.5, summation([a,i,b,j,c,k,d,l]), lcc0, rcc0, tcc1, bra1, comm, kHF])

rcc1 = tensor([b,j,c,k], 'R^{(1)}')
comm = commutator(phi, el_oper(b, j) * el_oper(c, k))

term2 = term(elmts=[0.5, summation([a,i,b,j,c,k]), lcc0, rcc1, bra1, comm, kHF])

exp = term1 + term2

# third order correction one term fundamentally different than for CIS(D)
rcc0 = tensor([c,k], 'R^\\text{CCS}')
c1 = commutator(phi, el_oper(b, j) )
comm = commutator(c1, el_oper(c, k))
tcc2 = tensor([b,j], 't^{(2)}')

exp =  term(elmts=[summation([a,i,b,j,c,k]), lcc0, rcc0, tcc2, bra1, comm, kHF])

c1 = commutator(phi, el_oper(c, k) )
c2 = commutator(c1, el_oper(d, l))
comm = commutator(c2, el_oper(e, m))
exp =  term(elmts=[summation([c,k,d,l,e,m]), bra2, comm, kHF])

## RHS for second-order singles ground state
#
#tcc1 = tensor([b,j,c,k], 't^{(1)}')
#comm = commutator(phi, el_oper(b, j) * el_oper(c, k))
#exp = term(sign='-', elmts=[0.5, summation([b,j,c,k]), tcc1, bra1, comm, kHF])
#
## expression for second-order doubles ground state amplitudes
#tcc1 = tensor([c,k,d,l], 't^{(1)}')
#comm = commutator(phi, el_oper(c, k) * el_oper(d, l))
#exp = term( elmts=[0.5, summation([c,k,d,l]), tcc1, bra2, comm, kHF] )
#
##comm = commutator(phi, el_oper(d, l))
##exp = term( elmts=[summation([c,k,d,l]), tcc1, bra2, el_oper(c,k), comm, kHF] )
#
## expression for second-order doubles transition amplitudes (first term)
#rcc0 = tensor([e,m], 'R^\\text{CCS}')
#tcc1 = tensor([c,k,d,l], 't^{(1)}')
#c1 = commutator(phi, el_oper(d, l) )
#comm = commutator(c1, el_oper(e, m))
#exp = term( elmts=[summation([c,k,d,l,e,m]), tcc1, rcc0, bra2, el_oper(c,k), comm, kHF] )

simple = exp.simplify()

result = simple.evaluate()



# create output file
filename = "output.tex"
out = open(filename, 'w')

out.write("\documentclass{article}\n")
out.write("\usepackage{amsmath}\n")
out.write("\usepackage{braket}\n")
out.write("\\begin{document}\n")

out.write('$'+exp.output()+' = $ \\\\ ~ \\\\ \n')
out.write('$ = '+simple.output()+' = $ \\\\ ~ \\\\ \n')
out.write('$ = '+result.output()+'$ \\\\ ~ \\\\ \n')

out.write("\end{document}\n")

