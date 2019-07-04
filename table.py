import numpy as np
import sys
import csv
import os

if len(sys.argv) < 2:
    print 'This script requires two arguments'
    print '\t1. directory with .out files from s mitigator programs'
    print '\t2. decimal number for data rounding'
    sys.exit()

a = np.zeros((7,3), dtype=np.float64)
mi = np.zeros((7,3), dtype=np.float64)
ma = np.zeros((7,3), dtype=np.float64)
x = np.zeros((7,3), dtype=np.float64)
decimals = int(sys.argv[2])

files=[]
lines=[]
N=7
k=0
for f in sorted(os.listdir(sys.argv[1])):
    if f.endswith(".out"):
        files.append(f)
        print 'working on', f, 'file'
        with open(f) as csvfile:
            lines = csvfile.readlines()
        x = np.genfromtxt(lines[-N:],delimiter=',')
        x = np.absolute(x)
        for i in xrange(7):
            a[i,k] = x[i,0]
            mi[i,k]= x[i,1]
            ma[i,k]= x[i,2]
        k+=1

a = np.around(a, decimals)
mi = np.around(mi, decimals)
ma = np.around(ma, decimals)

print '\\begin{table}'
print '\t\\centering'
print '\t\\begin{tabular}{cccccc}'
print '\t\t& \multicolumn{2}{c}{$m>\mathtext{EW}$} &\
       \multicolumn{2}{c}{$\Delta m^2\gtrsim 100{\mathtext{eV}}^2$} &\
       \multicolumn{2}{c}{$\Delta m^2\sim 0.1-1{\mathtext{eV^2}}^2$}\\'
print '\t\t\midrule'
print '\t\t\tEntry & $T$ & $\mathcal{T}$ & $T$ & $\mathcal{T}$ & $T$ &\
       $\mathcal{T}$'
print '\t\t\midrule'
print '\t\t\t$(1,1)$ & $', a[0,0], '$ & $', mi[0,0], '\div', ma[0,0], '$ &\
                       $', a[0,1], '$ & $', mi[0,1], '\div', ma[0,1], '$ &\
                       $', a[0,2], '$ & $', mi[0,2], '\div', ma[0,2], '$ \\'
print '\t\t\t$(2,2)$ & $', a[2,0], '$ & $', mi[2,0], '\div', ma[2,0], '$ &\
                       $', a[2,1], '$ & $', mi[2,1], '\div', ma[2,1], '$ &\
                       $', a[2,2], '$ & $', mi[2,2], '\div', ma[2,2], '$ \\'
print '\t\t\t$(3,3)$ & $', a[5,0], '$ & $', mi[5,0], '\div', ma[5,0], '$ &\
                       $', a[5,1], '$ & $', mi[5,1], '\div', ma[5,1], '$ &\
                       $', a[5,2], '$ & $', mi[5,2], '\div', ma[5,2], '$ \\'
print '\t\t\t$(2,1)$ & $', a[1,0], '$ & $', mi[1,0], '\div', ma[1,0], '$ &\
                       $', a[1,1], '$ & $', mi[1,1], '\div', ma[1,1], '$ &\
                       $', a[1,2], '$ & $', mi[1,2], '\div', ma[1,2], '$ \\'
print '\t\t\t$(3,1)$ & $', a[2,0], '$ & $', mi[2,0], '\div', ma[2,0], '$ &\
                       $', a[2,1], '$ & $', mi[2,1], '\div', ma[2,1], '$ &\
                       $', a[2,2], '$ & $', mi[2,2], '\div', ma[2,2], '$ \\'
print '\t\t\t$(3,2)$ & $', a[3,0], '$ & $', mi[3,0], '\div', ma[3,0], '$ &\
                       $', a[3,1], '$ & $', mi[3,1], '\div', ma[3,1], '$ &\
                       $', a[3,2], '$ & $', mi[3,2], '\div', ma[3,2], '$ \\'
print '\t\t\midrule'
print '\t\t\t$\sigma_3$ & $', ma[6,0], '$ & $', ma[6,1], '$ & $', ma[6,1],\
       '$\\'
print '\t\t\\bottomrule'
print '\t\end{tabular}'
print '\end{table}'
