import random
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import logging
import collections
DEVICE = 'cpu'

#Set up logging format
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s: %(message)s')

class DQNAgent(torch.nn.Module):
    def __init__(self, params):
        super().__init__()
        self.reward = 0
        self.gamma = 0.9
        self.short_memory = np.array([])
        self.agent_target = 1
        self.agent_predict = 0
        self.learning_rate = params['learning_rate']        
        self.epsilon = 1
        self.actual = []
        self.first_layer = params['first_layer_size']
        self.second_layer = params['second_layer_size']
        self.third_layer = params['third_layer_size']
        self.memory = collections.deque(maxlen=params['memory_size'])
        self.weights = params['weights_path']
        self.load_weights = params['load_weights']
        self.optimizer = None
        self.network()

    def network(self):
        # Layers
        self.f1 = nn.Linear(4, self.first_layer)
        self.f2 = nn.Linear(self.first_layer, self.second_layer)
        self.f3 = nn.Linear(self.second_layer, self.third_layer)
        self.f4 = nn.Linear(self.third_layer, 8)
        # weights
        if self.load_weights:
            self.model = self.load_state_dict(torch.load(self.weights))
            print("weights loaded")

    def forward(self, x):
        x = F.relu(self.f1(x))
        x = F.relu(self.f2(x))
        x = F.relu(self.f3(x))
        x = F.softmax(self.f4(x), dim=-1)
        return x

    def get_state(self, environment):
        """
        Return the state which is a numpy array of (array dimension) values representing:
            - +(-)1 if the element 11 is above(below) limit, 0 otherwise
            - +(-)1 if the element 12 is above(below) limit, 0 otherwise
            .
            .
            .
        """
        state = []
        state_index = 0
        for row in range(environment.limit_matrix_dimensions[0]):
            for column in range(environment.limit_matrix_dimensions[1]):
                if environment.result_matrix[row,column] < environment.limit_matrix[row,column]:
                    state.append(-1)
                elif environment.result_matrix[row,column] > environment.limit_matrix[row,column]:
                    state.append(+1)
                else:
                    state.append(0)
                state_index += 1
        logging.debug("State list %s", state)
        return np.asarray(state)

    def set_reward(self, result_matrix, limit_matrix): 
        self.reward = 0
        record = np.sum(np.square(np.subtract(result_matrix, limit_matrix)))
        if record < 0.5:
            self.reward = 10
        else:
            self.reward = 0
        return self.reward
