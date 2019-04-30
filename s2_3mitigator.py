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
expmatrix1 = [[1.3e-3, 0, 0], [6.8e-4, 2.2e-4, 0], [2.7e-3, 1.2e-3, 2.8e-3]]
expmatrix2 = [[2.4e-2, 0, 0], [2.5e-2, 2.2e-2, 0], [6.9e-2, 1.2e-2, 1.0e-1]]
expmatrix3 = [[1.0e-2, 0, 0], [1.7e-2, 1.4e-2, 0], [4.5e-2, 5.3e-2, 1.0e-1]]

exp_boundry = int(sys.argv[2])

if exp_boundry == 1:
    compmatrix = np.eye(3) - expmatrix1
    s = [1.0, 0.9999, 0.9999]
    sjump = 0.0001
    decraesequantity = 10000
elif exp_boundry == 2:
    compmatrix = np.eye(3) - expmatrix2
    s = [1.0, 0.999, 0.999]
    sjump = 0.001
    decraesequantity = 1000
elif exp_boundry == 3:
    compmatrix = np.eye(3) - expmatrix3
    s = [1.0, 0.999, 0.999]
    sjump = 0.001
    decraesequantity = 1000
else:
    decraesequantity = 1
    compmatrix = np.eye(3) - expmatrix3
    s = [1.0, 1.0, 1.0]
    sjump = 0.001

loopbreak = 1
control = 1
nthrows = int(sys.argv[3])

for i in range(decraesequantity):
    maj_problem = 0
    positive_fit = 0
    negative_fit = 0
    for j in tqdm(range(nthrows)):
        e = [uniform(compmatrix[0,0], 1), uniform(compmatrix[1,1], 1),
             uniform(compmatrix[2,2], 1)]
        result = create_alpha(s, e, 0)
        if result.all() == 1:
            maj_problem += 1
        else:
            if (compmatrix[0,0] <= result[0,0] <= 1 and
                compmatrix[1,1] <= result[1,1] <= 1 and
                compmatrix[2,2] <= result[2,2] <= 1 and
                0 < np.abs(result[1,0]) <= np.abs(compmatrix[1,0]) and
                0 < np.abs(result[2,0]) <= np.abs(compmatrix[2,0]) and
                0 < np.abs(result[2,1]) <= np.abs(compmatrix[2,1])):
                positive_fit += 1
                f.write('----------------------------------------\n')
                string = 'singular values: ' + str(s) + '\neigen values: ' +\
                str(e) + '\n'
                f.write(string)
                f.write('created matrix:\n')
                f.write(str(result) + '\n')
                eigv, vec = LA.eig(result)
                uh, singv ,vh = LA.svd(result)
                f.write('singular values of produced matrix: ' +\
                        str(np.around(singv, decimals=3)) + '\n')
                f.write('eigenvalues of produced matrix: ' + str(eigv) + '\n')
            else:
                negative_fit += 1
        if positive_fit == loopbreak:
            break
    f.write('==========')
    f.write(' Majorization problem: ' + str(maj_problem) + '/' + str(nthrows)
            + ' for s2: ' + str(s[1]) + ' and s3: ' + str(s[2]) + '\n')
    if positive_fit == 0:
        f.write('Could not find any matrix for given criteria.\n')
        break
    elif s[1] == sjump and s[2] == sjump:
        f.write('Achieved lowest singular value: ' + str(s[2]) + ' break\n')
        break
    elif control == 1:
        s[2] -= sjump
        control = 2
    elif control == 2:
        s[1] -= sjump
        control = 1
f.close()

