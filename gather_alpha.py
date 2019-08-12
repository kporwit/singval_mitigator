import numpy as np
import csv
from alpha_utils import create_alpha
import numpy.linalg as LA


expmatrix1 = np.array([[1.3e-3, 0, 0],\
                       [6.8e-4, 2.2e-4, 0],\
                       [2.7e-3, 1.2e-3, 2.8e-3]], dtype=np.float64)
expmatrix2 = np.array([[2.4e-2, 0, 0],\
                       [2.5e-2, 2.2e-2, 0],\
                       [6.9e-2, 1.2e-2, 1.0e-1]], dtype=np.float64)
expmatrix3 = np.array([[1.0e-2, 0, 0],\
                       [1.7e-2, 1.4e-2, 0],\
                       [4.5e-2, 5.3e-2, 1.0e-1]], dtype=np.float64)


compmatrix = np.eye(3) - expmatrix3

eigen_values = np.zeros(3,dtype=np.float64)
singular_values = np.array([1.00,1.00,0.899],dtype=np.float64)
positive_fit, maj_problem, negative_fit = 0, 0, 0
decimals=5

result_filename = 'gather_alpha_1.csv'
output = open(result_filename, 'w')
output.write('''alpha[1,1], alpha[2,2], alpha[3,3], alpha[2,1], alpha[3,1], \
alpha[3,2], s[1], s[2], s[3]\n''')

with open('major_1.csv') as csvfile:
    csv = csv.reader(csvfile, delimiter=',')
    for row in csv:
        eigen_values = row
        eigen_values = map(np.float64, eigen_values)
        created_alpha = create_alpha(singular_values,
                                     eigen_values,
                                     0)
        if (created_alpha==1).all():
            maj_problem += 1
        else:
            eigv, vec = LA.eig(created_alpha)
            u, singv ,vh = LA.svd(created_alpha)
            if (sorted(np.around(singv, decimals)) == sorted(singular_values) and
                       compmatrix[0,0] <= created_alpha[0,0] <= 1 and
                       compmatrix[1,1] <= created_alpha[1,1] <= 1 and
                       compmatrix[2,2] <= created_alpha[2,2] <= 1 and
                       np.abs(created_alpha[1,0]) <= np.abs(compmatrix[1,0]) and
                       np.abs(created_alpha[2,0]) <= np.abs(compmatrix[2,0]) and
                       np.abs(created_alpha[2,1]) <= np.abs(compmatrix[2,1])):
                positive_fit += 1
                output_string = '%f, %f, %f, %f, %f, %f, %f, %f, %f'\
                                %(created_alpha[0,0], created_alpha[1,1],\
                                  created_alpha[2,2], created_alpha[1,0],\
                                  created_alpha[2,0], created_alpha[2,1],\
                                  singv[0], singv[1], singv[2])
                output.write(output_string + '\n')
            else:
                negative_fit += 1
print 'majorization problems:\n', maj_problem
print 'matrices within bounds:\n', positive_fit
print 'matrices outside bounds:\n', negative_fit
output.close()


