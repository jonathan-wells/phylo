#!/usr/bin/env python3

import numpy as np
import time

class HMM(object):
    """Base class for HMMs."""

    def __init__(self, nstates, initprobs, transition, emission):
        self.states = np.array([i for i in range(nstates)])
        self.initprobs = initprobs
        self.transition = transition  # Transition probs between hidden states
        self.emission = emission  # Emission probs given hidden states
        self.emission_states = np.array([i for i in range(emission.shape[1])])
        self._check_input()
    
    def _check_input(self, err=1e-5):
        if self.initprobs.shape != self.states.shape:
            raise ValueError('incorrect initprob format')
        if self.initprobs.shape[0] != self.transition.shape[0]:
            raise ValueError()
        if self.emission.shape[0] != self.transition.shape[0]:
            raise ValueError()
        if (1 - sum(self.initprobs)) > err:
            raise ValueError('initprobs do not sum to 1')
        for i in self.states:
            if abs(1 - sum(self.transition[i])) > err:
                raise ValueError('transition probs do not sum to 1')
        for i in self.states:
            if abs(1 - sum(self.emission[i])) > err:
                raise ValueError('emission probs do not sum to 1')


    def simulate(self, t):
        simhidden = [np.random.choice(self.states, p=self.initprobs)]
        for i in range(1, t):
            pred = np.random.choice(self.states, 
                                    p=self.transition[simhidden[i-1], :])
            simhidden.append(pred)
        
        simobs = []
        for i in simhidden:
            pred = np.random.choice(self.emission_states,
                                    p=self.emission[i, :])
            simobs.append(pred)
        return simobs, simhidden

    def viterbi(self, observations):
        log_initprobs = np.log(self.initprobs)
        log_transition = np.log(self.transition)
        log_emission = np.log(self.emission)
    
        tracemat = np.zeros((self.states.shape[0], len(observations)), int)
        scoremat = np.zeros((self.states.shape[0], len(observations)))
        
        # Populate matrices
        for i in self.states:
            scoremat[i, 0] = log_initprobs[i] + \
                             log_emission[i, observations[0]]
        
        def prob(i, j, k, obs):
            return scoremat[k, j-1] + \
                   log_transition[k, i] + \
                   log_emission[i, obs]

        for j, obs in enumerate(observations[1:], 1):
            for i in self.states:
                probs = np.array([prob(i, j, k, obs) for k in self.states])
                scoremat[i, j] = probs.max()
                tracemat[i, j] = np.argmax(probs)
         
        # Traceback 
        z = np.argmax(scoremat[:, -1])
        preds = []
        scores = []
        j = len(observations) - 1
        while j >= 0:
            preds.append(z)
            scores.append(scoremat[z, j])
            z = tracemat[z, j]
            j -= 1
        
        return preds[::-1], scores[::-1]
                

def simulate_and_test(n):
    state_code = ['G', 'I']
    obs_code = ['A', 'C', 'T', 'G']
    nstates = len(state_code)
    initprobs = np.array([0.99, 0.01], float)
    transition = np.array([[0.99, 0.01],
                           [0.05, 0.95]], float)
    emission = np.array([[0.4, 0.1, 0.4, 0.1],
                         [0.2, 0.3, 0.2, 0.3]], float)
    test = HMM(nstates, initprobs, transition, emission)
    so, sh = test.simulate(n)
    preds, scores = test.viterbi(so)
    matches = ''
    for i, val in enumerate(sh):
        if val == preds[i]:
            matches += '|'
        else:
            matches += ' '
    print('observations:', ''.join([obs_code[i] for i in so]))
    print()
    print('states:      ', ''.join([state_code[i] for i in sh]))
    print('             ', matches),
    print('predictions: ', ''.join([state_code[i] for i in preds]))
    print(scores)

if __name__ == '__main__':
    simulate_and_test(90)
