

```julia
using QuantumCircuits
using QuantumCircuits.QML
using QuantumCircuits.QCircuits.Qiskit
using QuantumCircuits.QCircuits.Qiskit: qiskit
using QuantumCircuits.QCircuits.Circuit
using QuantumCircuits.QCircuits.Circuit: toQiskit, getCode
using QuantumCircuits.QCircuits.Gates: CX
using QuantumCircuits.Execute: generate_mesuere_circuits, extractProbability, correctMeasures
```

# 2 Introduction

Our goal is to implement the Trotterization to evolve the specified state $$|110\rangle$$, under the specified Hamiltonian, for the specified duration $$\pi$$ but utilize as many Trotter steps as is possible. The algorithm, in principle, should be executed for any state, Hamiltonian, and duration.

Unfortunately on the current quantum computer, this is impossible regarding the noises. Therefore I've proposed an algorithm that can break the Trotterization into pieces with a shorter depth which can be executed on current quantum computers. The main points and requirements of the algorithm:
- It breaks the Trotterization into pieces with a shorter depth.
- It is convergent to the final state.
- The designed algorithm should allow us to find the simulated state for an arbitrary number of qubits.

The first stage is to find a way how to effectively implement a trotter step, the details of that can be found in the notebook [03_trotter_step](03_trotter_step.ipynb). 

# 3 The overview of algorithm

Let's assume, that we would like to implement an algorithm with 10 Trotter steps, but using at most 2 Trotter steps in single circuit execution. We can do this by implementing the circuit below and using the gradient descent method to find the parameters $\theta_1$ that minimize our cost function. The cost function is chosen that after finding the optimal parameters the output status is $|0\rangle_3$.
Why I use $|000\rangle$ as a target state, there are a few reasons:
1. We need to avoid any superposition and use only base state $|000\rangle$, $|001\rangle$, eg. Using superposition implies variance > 0 during expectation measurement even on the perfect device, for the base state the variance is equal to 0. This is important because any fluctuation impact derivative calculation using the shift-parameter rule.
1. $|000\rangle$ state in natural choice because this is the starting point for quantum device.

First, we would like to fit the expected distribution. The use of entropy loos is a natural choice point. In our case, the expected distribution has only one non zero value for state 000, so the entropy loss looks like: $-\ln{P(|000\rangle)}$ But the other requirement is that the error for the other state should be equally distributed, so we prefer $\frac{|001\rangle + |010\rangle}{\sqrt{2}}$ to $|001\rangle$ but the entropy for both of the states is equal. Therefore our loss function is equal to: $-\ln{P(|000\rangle)} + \sum_{i = |001\rangle}^{|111\rangle}P(i)^2$ (see [04_algorithm_evaluation](04_algorithm_evaluation.ipynb))
