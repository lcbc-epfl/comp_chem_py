#!/usr/bin/env python
from subprocess import check_call

# create define input file
lines = """

a coord
*
no
b all aug-cc-pVDZ
*
eht
y
1
y
scf
iter
100

dft
on
func slater-dirac-exchange
*
ex
ciss
*
a 5
*
*
y
*
""".split('\n')

with open('define.inp', 'w') as df:
    for line in lines:
        df.write(line+'\n')


# execute define to get control file
check_call(['bash','-c', 'define < define.inp'])

# get running state from file SH_STATE.dat (Only when surface hopping)
try:
    with open('../SH_STATE.dat', 'r') as state:
        last =  state.read().split('\n')[-2]
        current = int(last.split()[-1])
        print('RUN ON SURFACE {}'.format(current))

except(IOError):
    current = 2 # state for which the gradient is calculated

# modify control file
with open('control', 'r') as ctrl:
    lines = []
    for line in ctrl:
        if '$soes' in line:
            lines.append('$exopt {}\n'.format( current ))
        lines.append(line)

with open('control', 'w') as ctrl:
    for line in lines:
        ctrl.write(line)

