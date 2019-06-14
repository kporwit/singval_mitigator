from alpha_utils import create_alpha
from random import seed, randrange, uniform
from tqdm import tqdm
from numpy import linalg as LA
import numpy as np
import time
import sys
import os.path

seed(time.time())

if len(sys.argv) < 3:
    print 'This script requires 3 arguments:'
    print '\t1: name of an output file,'
    print '\t2: number of experimental matrix,'
    print '\t3: number of random throws.'
    sys.exit()

print 'Running s123mitigator...'

filename = str(sys.argv[1])
if os.path.isfile(filename) == 1:
    print 'File exist! Exit.'
    sys.exit()
else:
    f = open(filename, 'w')

f.write('========================================\n' +
        str(time.strftime('%d/%m/%Y|%H:%M:%S')))
f.write('\n========================================\n' + str(sys.argv) +
        '\n========================================\n')

#comparision betwen experimental boundries
expmatrix1 = np.array([[1.3e-3, 0, 0],\
                       [6.8e-4, 2.2e-4, 0],\
                       [2.7e-3, 1.2e-3, 2.8e-3]], dtype=np.float64)
expmatrix2 = np.array([[2.4e-2, 0, 0],\
                       [2.5e-2, 2.2e-2, 0],\
                       [6.9e-2, 1.2e-2, 1.0e-1]], dtype=np.float64)
expmatrix3 = np.array([[1.0e-2, 0, 0],\
                       [1.7e-2, 1.4e-2, 0],\
                       [4.5e-2, 5.3e-2, 1.0e-1]], dtype=np.float64)

#matrix for proper final result
fin_result = np.zeros((3,3), dtype=np.float64)
#matrices for maximal and minimal entries of generated alpha
maximal = np.zeros((3,3), dtype=np.float64)
minimal = np.zeros((3,3), dtype=np.float64)
#variable for latex output
tex_out = 0

exp_boundry = int(sys.argv[2])

if exp_boundry == 1:
    compmatrix = np.eye(3) - expmatrix1
    s = np.array([1.0, 0.9999, 0.9999], dtype=np.float64)
    sjump = np.double(0.0001)
    decraesequantity = 10000
    decimals = 4
elif exp_boundry == 2:
    compmatrix = np.eye(3) - expmatrix2
    s = np.array([1.0, 0.999, 0.999], dtype=np.float64)
    sjump = np.double(0.001)
    decraesequantity = 1000
    decimals = 3
elif exp_boundry == 3:
    compmatrix = np.eye(3) - expmatrix3
    s = np.array([1.0, 0.999, 0.999], dtype=np.float64)
    sjump = np.double(0.001)
    decraesequantity = 1000
    decimals = 3

loopbreak = 1
control = 1
nthrows = int(sys.argv[3])

for i in range(decraesequantity):
    maj_problem = 0
    positive_fit = 0
    negative_fit = 0
    #for j in tqdm(range(nthrows)):
    for j in range(nthrows):
        e = np.array([uniform(compmatrix[0,0], 1), uniform(compmatrix[1,1], 1),
             uniform(compmatrix[2,2], 1)], dtype=np.float64)
        result = create_alpha(s, e, 0)
        if result.all() == 1:
            maj_problem += 1
        else:
            eigv, vec = LA.eig(result)
            uh, singv ,vh = LA.svd(result)
            if (sorted(np.around(singv, decimals)) == sorted(s) and
                compmatrix[0,0] <= result[0,0] <= 1 and
                compmatrix[1,1] <= result[1,1] <= 1 and
                compmatrix[2,2] <= result[2,2] <= 1 and
                0 < np.abs(result[1,0]) <= np.abs(compmatrix[1,0]) and
                0 < np.abs(result[2,0]) <= np.abs(compmatrix[2,0]) and
                0 < np.abs(result[2,1]) <= np.abs(compmatrix[2,1])):
                positive_fit += 1
                tex_out = 1
                absresult = np.absolute(result)
                if i == 0:
                    maximal = np.copy(absresult)
                    minimal = np.copy(absresult)
                np.putmask(maximal, absresult>maximal, absresult)
                np.putmask(minimal, absresult<minimal, absresult)
                f.write('----------------------------------------\n')
                string = 'singular values: ' + str(s) + '\neigen values: ' +\
                str(e) + '\n'
                f.write(string)
                f.write('created matrix:\n')
                f.write(str(result) + '\n')
                f.write('rounded singular values of produced matrix: ' +\
                        str(np.around(singv, decimals)) + '\n')
                f.write('singular values of produced matrix: ' +\
                        str(singv) + '\n')
                f.write('eigenvalues of produced matrix: ' + str(eigv) + '\n')
                fin_result = np.copy(result)
            else:
                negative_fit += 1
        if positive_fit == loopbreak:
            break
    f.write('==========')
    f.write(' Majorization problem: ' + str(maj_problem) + '/' + str(nthrows)
            + 'for s1: ' + str(s[0]) + ' for s2: ' + str(s[1])
            + ' and s3: ' + str(s[2]) + '\n')
    if positive_fit == 0:
        f.write('Could not find any matrix for given criteria.\n')
        break
    elif s[2] == sjump:
        f.write('Achieved lowest singular value: ' + str(s[2]) + ' break\n')
        break
    elif control == 1:
        s[2] -= sjump
        s[2] = np.around(s[2], decimals)
        control = 2
    elif control == 2:
        s[1] -= sjump
        s[1] = np.around(s[1], decimals)
        control = 3
    elif control ==3:
        s[0] -=sjump
        s[0] = np.around(s[0], decimals)
        control = 1
if tex_out == 1:
    np.around(fin_result, decimals=5)
    np.around(eigv, decimals=5)
    np.around(singv, decimals=3)
    f.write('----------------------------------------\nMaximal values:\n' + str(maximal) + '\nMinimal\
            values:\n' + str(minimal) + '\n')

    f.write('-----Latex script with output-----\n\n')
    if exp_boundry == 1:
        f.write('$$T_{\unit{m>246}{GeV}}=\n')
    elif exp_boundry == 2:
        f.write('$$T_{{\unit{\Delta m^2\gtrsim100}{eV^2}}}=\n')
    else:
        f.write('$$T_{\unit{\Delta m^2 \sim0.1-1}{eV^2}}=\n')
    f.write(r'\begin{pmatrix}')
    f.write('\n\t' + str(fin_result[0,0]) + ' & ' + str(fin_result[0,1])\
            + ' & ' + str(fin_result[0,2]) + ' \\\\\n')
    f.write('\t' + str(fin_result[1,0]) + ' & ' + str(fin_result[1,1])\
            + ' & ' + str(fin_result[1,2]) + ' \\\\\n')
    f.write('\t' + str(fin_result[2,0]) + ' & ' + str(fin_result[2,1])\
            + ' & ' + str(fin_result[2,2]) + ' \\\\\n\end{pmatrix}$$\n')
    if exp_boundry == 1:
        f.write('$$\sigma(T_{\unit{m>246}{GeV}}) = \{' +\
                str(np.around(singv[0], decimals)) + ', ' +\
                str(np.around(singv[1], decimals)) + ', ' +\
                str(np.around(singv[2], decimals)) +\
                '\}$$\n')
        f.write('$$\lambda(T_{\unit{m>246}{GeV}}) = \{' +\
                str(eigv[0]) + ', ' + str(eigv[1]) + ', ' + str(eigv[2]) +\
                '\}$$\n')
    elif exp_boundry == 2:
        f.write('$$\sigma(T_{\unit{\Delta m^2\gtrsim100}{eV^2}}) = \{' +\
                str(np.around(singv[0], decimals)) + ', ' +\
                str(np.around(singv[1], decimals)) + ', ' +\
                str(np.around(singv[2], decimals)) +\
                '\}$$\n')
        f.write('$$\lambda(T_{\unit{\Delta m^2\gtrsim100}{eV^2}}) = \{' +\
                str(eigv[0]) + ', ' + str(eigv[1]) + ', ' + str(eigv[2]) +\
                '\}$$\n')
    else:
        f.write('$$\sigma(T_{\unit{\Delta m^2 \sim0.1-1}{eV^2}}) = \{' +\
                str(np.around(singv[0], decimals)) + ', ' +\
                str(np.around(singv[1], decimals)) + ', ' +\
                str(np.around(singv[2], decimals)) +\
                '\}$$\n')
        f.write('$$\lambda(T_{\unit{\Delta m^2 \sim0.1-1}{eV^2}}) = \{' +\
                str(eigv[0]) + ', ' + str(eigv[1]) + ', ' + str(eigv[2]) +\
                '\}$$\n')
    f.write(r'$$\begin{pmatrix}')
    f.write('\n\t' + str(minimal[0,0]) + ' - ' + str(maximal[0,0]) + ' & 0 & 0 \\\\\n')
    f.write('\t' + str(minimal[1,0]) + ' - ' + str(maximal[1,0]) + ' & '\
                 + str(minimal[1,1]) + ' - ' + str(maximal[1,1]) + ' & 0 \\\\\n')
    f.write('\t' + str(minimal[2,0]) + ' - ' + str(maximal[2,0]) + ' & '\
                 + str(minimal[2,1]) + ' - ' + str(maximal[2,1]) + ' & '\
                 + str(minimal[2,2]) + ' - ' + str(maximal[2,2]) + '\\\\\n')
    f.write('\end{pmatrix}$$\n')
f.close()
