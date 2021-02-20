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

        logging.info("Provided singular values %s and limit matrix \n %s", self.singular_values, self.limit_matrix)
        self.singular_matrix = np.diag(self.singular_values)
        logging.debug("Singular matrix:\n %s", self.singular_matrix)

    def generate_unitary_matrix(self, dimension):
        #for now this matrix is not unitary
        unitary_matrix = np.random.rand(dimension, dimension)
        return unitary_matrix

    def generate_result_matrix(self):
        U = self.generate_unitary_matrix(self.singular_values_dimension)
        logging.debug("Generated U matrix \n%s", U)
        V = self.generate_unitary_matrix(self.singular_values_dimension)
        logging.debug("Generated V matrix \n%s", V)
        result_matrix = np.matmul(np.matmul(U, self.singular_matrix), V)
        logging.debug("Generated result matrix \n%s", result_matrix)
        return result_matrix

singular_values = np.array([1, 0.99], dtype=float)
limit_matrix = np.array([[1,2],[1,2]], dtype=float)
test = Environment(singular_values, limit_matrix)
test.generate_result_matrix()
