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
        self.reward = 1000
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
                state.append(environment.limit_matrix[row,column] - environment.result_matrix[row,column])
                #if environment.result_matrix[row,column] < environment.limit_matrix[row,column]:
                #elif environment.result_matrix[row,column] > environment.limit_matrix[row,column]:
                #    state.append(+1)
                #else:
                #    state.append(0)
                state_index += 1
        logging.debug("State list %s", state)
        return np.asarray(state)

    def replay_new(self, memory, batch_size):
        """
        Replay memory
        """
        if len(memory) > batch_size:
            minibatch = random.sample(memory, batch_size)
        else:
            minibatch = memory
        for state, u_action, v_action, reward, next_state, done in minibatch:
            self.train()
            torch.set_grad_enabled(True)
            target = reward
            next_state_tensor = torch.tensor(np.expand_dims(next_state, 0), dtype=torch.float32).to(DEVICE)
            state_tensor = torch.tensor(np.expand_dims(state, 0), dtype=torch.float32, requires_grad=True).to(DEVICE)
            if not done:
                target = reward + self.gamma * torch.max(self.forward(next_state_tensor[0]))
            output = self.forward(state_tensor)
            target_f = output.clone()
            target_f[0][np.argmax(np.vstack((u_action, v_action)))] = target
            target_f.detach()
            self.optimizer.zero_grad()
            loss = F.mse_loss(output, target_f)
            loss.backward()
            self.optimizer.step()

    def remember(self, state, u_action, v_action, reward, next_state, done):
        """
        Store the <state, action, reward, next_state, is_done> tuple
        in a memory buffer for replay memory
        """
        self.memory.append((state, u_action, v_action, reward, next_state, done))

    def set_reward(self, result_matrix, limit_matrix): 
        old_reward = self.reward
        reward = np.sum(np.square(np.subtract(result_matrix, limit_matrix)))
        if reward < old_reward:
            self.reward = reward
        else:
            self.reward = -1
        logging.debug("Record is %s which corresponds to reward %s", reward, self.reward)
        return self.reward

    def train_short_memory(self, state, u_action, v_action, reward, next_state, done):
        """
        Train the DQN agent on the <state, action, reward, next_state, is_done>
        tuple at the current timestep.
        """
        self.train()
        torch.set_grad_enabled(True)
        target = reward
        next_state_tensor = torch.tensor(next_state.reshape((1,4)), dtype=torch.float32).to(DEVICE)
        state_tensor = torch.tensor(state.reshape((1,4)), dtype=torch.float32, requires_grad=True).to(DEVICE)
        if not done:
            target = reward + self.gamma * torch.max(self.forward(next_state_tensor[0]))
        #logging.debug("Short memory params: target %s next_state_tensor %s state_tensor %s", target, next_state_tensor, state_tensor)
        output = self.forward(state_tensor)
        target_f = output.clone()
        target_f[0][np.argmax(np.vstack((u_action, v_action)))] = target
        target_f.detach()
        self.optimizer.zero_grad()
        loss = F.mse_loss(output, target_f)
        #logging.debug("Short memory NN params: output %s target_f %s loss %s ", output, target_f, loss)
        loss.backward()
        self.optimizer.step()
