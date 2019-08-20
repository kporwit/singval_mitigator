import numpy as np
import csv
import sys
import numpy.linalg as LA


def create_lt_matrix(lt_elements_list):
    lte = lt_elements_list 
    lt_matrix = np.array([[lte[0], 0.0 ,0.0],\
                          [lte[3], lte[1], 0.0],\
                          [lte[4], lte[5], lte[2]]])
    return lt_matrix


input_filename = str(sys.argv[1])
minimal_alpha_values = np.ones(6)
maximal_alpha_values = np.full(6,-40.0)
minimal_unitary_mat = np.ones((3,3))
maximal_unitary_mat = np.zeros((3,3))
abs_minimal_unitary_mat = np.ones((3,3))
abs_maximal_unitary_mat = np.zeros((3,3))

with open(input_filename) as csvfile:
    csv = csv.reader(csvfile, delimiter=',')
    for row in csv:
        row = map(np.float, row)
        alpha_elements = row[0:6]
        singular_values = row[6:]
        input_matrix = create_lt_matrix(alpha_elements)
        minimal_alpha_values = np.where(minimal_alpha_values<alpha_elements,
                                        minimal_alpha_values, alpha_elements)
        maximal_alpha_values = np.where(maximal_alpha_values>alpha_elements,
                                        maximal_alpha_values, alpha_elements)
        minimal_alpha_matrix = create_lt_matrix(minimal_alpha_values)
        unitary_matrix, sv, vh = LA.svd(input_matrix)
        #unitary_matrix = np.around(unitary_matrix, decimals=6)
        minimal_unitary_mat = np.where(minimal_unitary_mat<unitary_matrix,
                                        minimal_unitary_mat, unitary_matrix)
        maximal_unitary_mat = np.where(maximal_unitary_mat>unitary_matrix,
                                        maximal_unitary_mat, unitary_matrix)
        abs_unitary_matrix = np.absolute(unitary_matrix)
        abs_minimal_unitary_mat = np.where(abs_minimal_unitary_mat<abs_unitary_matrix,
                                           abs_minimal_unitary_mat, abs_unitary_matrix)
        abs_maximal_unitary_mat = np.where(abs_maximal_unitary_mat>abs_unitary_matrix,
                                           abs_maximal_unitary_mat, abs_unitary_matrix)
        print 'gather_alpha output:\n', create_lt_matrix(row)
        print 'minimal values:\n', create_lt_matrix(minimal_alpha_values)
        print 'maximal values:\n', create_lt_matrix(maximal_alpha_values)
        print 'unitary matrix\n', unitary_matrix
        print 'minimal unitary:\n', minimal_unitary_mat
        print 'maximal unitary:\n', maximal_unitary_mat
        print 'absolute minimal unitary:\n', abs_minimal_unitary_mat
        print 'absolute maximal unitary:\n', abs_maximal_unitary_mat
        print '='*50


