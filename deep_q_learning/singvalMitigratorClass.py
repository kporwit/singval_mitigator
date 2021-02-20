import numpy as np
import logging
import sys

#Set up logging format
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s: %(message)s')

class Environment:
    """ Initialize environment """

    def __init__(self, singular_values, limit_matrix):
        self.singular_values = singular_values
        self.limit_matrix = limit_matrix

        if type(self.singular_values) is not np.ndarray:
            logging.error('Provided singular values are not in a numpy array form')
            sys.exit(1)

        self.limit_matrix_dimensions = self.limit_matrix.shape
        if self.limit_matrix_dimensions[0] != self.limit_matrix_dimensions[1]:
            logging.error('Provided limit matrix is not square matrix')
            sys.exit(2)

        self.singular_values_dimension = len(self.singular_values)
        if self.singular_values_dimension != self.limit_matrix_dimensions[0]:
            logging.error('Provided limit matrix dimensions are not aligned with the singular list dimension')
            sys.exit(2)

        logging.debug("Provided singular values %s and limit matrix \n %s", self.singular_values, self.limit_matrix)
        self.singular_matrix = np.diag(self.singular_values)

        self.U = np.random.rand(self.singular_values_dimension, self.singular_values_dimension)
        self.V = np.random.rand(self.singular_values_dimension, self.singular_values_dimension)
        logging.debug("Generated U matrix \n%s", self.U)
        logging.debug("Generated V matrix \n%s", self.V)

        self.result_matrix = self.generate_result_matrix()

    def perform_u_change(self, change_u):
        self.U = np.multiply(self.U, change_u)
        logging.debug("Generated U matrix \n%s", self.U)
        return self.U

    def perform_v_change(self, change_v):
        self.V = np.multiply(self.V, change_v)
        logging.debug("Generated V matrix \n%s", self.V)
        return self.V

    def generate_result_matrix(self):
        self.result_matrix = np.matmul(np.matmul(self.U, self.singular_matrix), self.V)
        logging.debug("Generated result matrix \n%s", self.result_matrix)
        return self.result_matrix

singular_values = np.array([1, 0.99], dtype=float)
limit_matrix = np.array([[1,2],[1,2]], dtype=float)
change_u = np.array([[1,1], [1,-1]], dtype=float)
change_v = np.array([[1,-1], [1,1]], dtype=float)
test = Environment(singular_values, limit_matrix)
test.perform_u_change(change_u)
test.perform_v_change(change_v)
test.generate_result_matrix()
print(test.U)
print(test.V)
print(test.result_matrix)
