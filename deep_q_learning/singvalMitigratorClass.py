import numpy as np
import logging
import sys
import argparse
from DQN import DQNAgent
from random import randint
import random
import statistics
import torch.optim as optim
import torch
import datetime

DEVICE = 'cpu'

def define_parameters():
        params = dict()
        # Neural Network
        params['epsilon_decay_linear'] = 1/100
        params['learning_rate'] = 0.00013629
        params['first_layer_size'] = 200    # neurons in the first layer
        params['second_layer_size'] = 20   # neurons in the second layer
        params['third_layer_size'] = 50    # neurons in the third layer
        params['episodes'] = 250          
        params['memory_size'] = 2500
        params['batch_size'] = 1000
        # Settings
        params['weights_path'] = 'weights/weights.h5'
        params['train'] = True
        params["test"] = False
        params['verbosity'] = 'DEBUG'

        return params


class Environment:
    """ Initialize environment """

    def __init__(self, singular_values, limit_matrix):
        self.singular_values = singular_values
        self.limit_matrix = limit_matrix
        self.change_U = 0
        self.change_V = 0

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
        #saves change matrix to the class atribute
        self.change_U = change_u
        #performs item wise multiplication on a U matrix
        self.U = np.multiply(self.U, change_u)
        logging.debug("Generated U matrix \n%s", self.U)
        return self.U

    def perform_v_change(self, change_v):
        self.change_V = change_v
        self.V = np.multiply(self.V, change_v)
        logging.debug("Generated V matrix \n%s", self.V)
        return self.V

    def generate_result_matrix(self):
        self.result_matrix = np.matmul(np.matmul(self.U, self.singular_matrix), self.V)
        logging.debug("Generated result matrix \n%s", self.result_matrix)
        return self.result_matrix

def run(params):
    singular_values = np.array([1, 0.99], dtype=float)
    limit_matrix = np.array([[0.934,0.174],[0.103,0.548]], dtype=float)
    change_u = np.array([[0.1,0.1], [0.1,-0.1]], dtype=float)
    change_v = np.array([[0.1,-0.1], [0.1,0.1]], dtype=float)

    #Initialize classes
    agent = DQNAgent(params)
    agent = agent.to(DEVICE)
    agent.optimizer = optim.Adam(agent.parameters(), weight_decay = 0, lr=params['learning_rate'])
    counter_games = 0
    record = 0
    total_score = 0
    while counter_games < params['episodes']:
        
        mitigation = Environment(singular_values, limit_matrix)

        mitigation.perform_u_change(change_u)
        mitigation.perform_v_change(change_v)
        mitigation.generate_result_matrix()
        agent.get_state(mitigation)

        steps = 0 #steps since the last positive reward
        while steps < 100:
            if not params['train']:
                agent.epsilon = 0.01
            else:
                agent.epsilon = 1 - (counter_games * params['epsilon_decay_linear'])
            #get old state
            state_old = agent.get_state(mitigation)

            #preform action based on agent.epsilon or choose the action
            if random.uniform(0,1) < agent.epsilon:
                #change U and V randomly
                final_u_move = np.multiply(change_u, np.random.uniform(-0.1,0.1,(2,2)))
                final_v_move = np.multiply(change_v, np.random.uniform(-0.1,0.1,(2,2)))
            else:
                #predict action based on the old state
                with torch.no_grad():
                    state_old_tensor = torch.tensor(state_old.reshape((1,4)), dtype=torch.float32).to(DEVICE)
                    prediction = agent(state_old_tensor)
                    final_u_move = prediction.detach().cpu().numpy()[0][0:4].reshape((2,2))
                    final_v_move = prediction.detach().cpu().numpy()[0][4:].reshape((2,2))
            #perform new move and get new state
            mitigation.perform_u_change(final_u_move)
            mitigation.perform_v_change(final_v_move)
            mitigation.generate_result_matrix()
            state_new = agent.get_state(mitigation)

            #set reward for the new state
            reward = agent.set_reward(mitigation.result_matrix, mitigation.limit_matrix)
            if reward > 0:
                steps = 0

            steps += 1

        logging.info("Game: %s reward %s", counter_games, agent.reward)
        counter_games += 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    params = define_parameters()
    #Set up logging format
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s: %(message)s')
    if params['train']:
        logging.info("Traininig...")
        params['load_weights'] = False
        run(params)
    if params['test']:
        logging.info("Testing...")
        params['load_weights'] = True
        run(params)
