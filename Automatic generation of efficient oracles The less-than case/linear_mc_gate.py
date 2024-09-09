

import numpy as np
from qiskit import QuantumCircuit
from collections import namedtuple




def _gate_U_pow(Ugate:np.array, denominator_exponent:int, inv:bool):
    '''
    Builds the quantum gate equivalent to Ugate^{1/param} either inverse or not.

    Input:
        - Ugate (np.array 2x2): Unitary matrix for single qubit transformation
        - param (int): power of the root
        - inv (bool): whether it is the inverse or direct matrix.
    
    Output:
        - Qiskit gate
    '''

    # Calculate exponent of the desired gate
    exp = 1 / np.abs(denominator_exponent)

    # Eigenvalues and eigenvectors of the Ugate
    eigenvals, eigenvects = np.linalg.eig(Ugate)

    value1, value2 = eigenvals[0], eigenvals[1]

    # Separate vectors
    v1 = eigenvects[:, [0]] # Needs to be a column, hence [:, [0]] and not [:, 0]
    v2 = eigenvects [:, [1]]

    # Conjugate vectors
    v1conj = v1.conj()
    v2conj = v2.conj()

    gate = np.power(value1 + 0j, exp) * np.dot(v1, v1conj.T) \
         + np.power(value2 + 0j, exp) * np.dot(v2, v2conj.T)
    
    # gate = np.power(value1 + 0j, exp) * v1 @ v1conj.T \
    #      + np.power(value2 + 0j, exp) * v2 @ v2conj.T

    if inv:
        gate = np.linalg.inv(gate)
        name = "Udg^1/%d"%denominator_exponent
    else:
        name = "U^1/%d"%denominator_exponent

    gate_circuit = QuantumCircuit(1, name = name)# '$\sqrt[%d]{U}$'%denominator_exponent)
    gate_circuit.unitary(gate, 0)

    return gate_circuit



def _build_P_gates(gate:np.array, nqubits:int, circuit:QuantumCircuit, inv:bool):
    '''
    Construction of the P gates.

    Input:
        - gate (np.array): single qubit gate to multi-control.
        It has to be a unitary matrix.
        - nqubits (int): numbr of qubits (including target).
        - circuit (QuantumCircuit): circuit on which the P gates are added.
        - inv (bool): whether it is the direct or the inverse.
    '''
    
    if inv:
        controls = list(range(1, nqubits-1))
    else:
        controls = list(range(nqubits-2, 0, -1))

    for k in controls:
        exponent = 2**(nqubits-1-k)
        gateU = _gate_U_pow(Ugate=gate, denominator_exponent=exponent, inv=inv)
        circuit.compose(gateU.control(1), [k, nqubits-1], inplace=True)


def _build_Q_gates(nqubits:int, circuit:QuantumCircuit, inv:bool):
    '''
    Construction of the Q gates.

    Input:
        - nqubits (int): numbr of qubits (including target).
        - circuit (QuantumCircuit): circuit on which the Q gates are added.
        - inv (bool): whether it is the direct or the inverse.
    '''

    qubit_pairs = namedtuple('qpairs', ['control', 'target'])

    if inv:
        start = 1
    else:
        start = 0

    qubit_pairs = [
        qubit_pairs(control, target) for target in range(nqubits) for control in range(start, target)
    ]

    qubit_pairs.sort(key=lambda k : k.control + k.target, reverse=True)

    for pair in qubit_pairs:
        exponent = pair.target - pair.control

        if pair.control == 0:
            exponent = exponent - 1
        
        param_turn = 2**exponent
        
        circuit.crx(theta=np.pi/param_turn, control_qubit=pair.control, target_qubit=pair.target)

    
    qubit_pairs = namedtuple('qpairs', ['control', 'target'])

    if inv:
        start = 0
    else:
        start = 1

    qubit_pairs = [
        qubit_pairs(control, target) for target in range(nqubits) for control in range(start, target)
    ]

    qubit_pairs.sort(key=lambda k : k.control + k.target, reverse=False)

    for pair in qubit_pairs:
        exponent = pair.target - pair.control

        if pair.control == 0:
            exponent = exponent - 1
        
        param_turn = 2**exponent

        circuit.crx(theta=-np.pi/param_turn, control_qubit=pair.control, target_qubit=pair.target)


def linear_mc_gate(gate:np.array, circuit:QuantumCircuit, controls:list, target: int):
    '''
    Construction of the linear multi-controlled gate.

    Input:
        - gate (np.array): single qubit gate to multi-control.
        It has to be a unitary matrix.
        - circuit (QuantumCircuit): circuit on which the gate is added.
        - controls (list): control qubits.
        - target (int): target qubit.
    '''

    nqubits = len(controls) + 1

    gate_circuit = QuantumCircuit(nqubits, name='MC_T%d'%target)

    _build_P_gates(gate=gate, nqubits=nqubits, circuit=gate_circuit, inv=False)

    exponent = 2**(nqubits-2)
    # print(exponent)
    gateU = _gate_U_pow(Ugate=gate, denominator_exponent=exponent, inv=False)
    gate_circuit.compose(gateU.control(1), [0, target], inplace=True)

    _build_Q_gates(nqubits=nqubits-1, circuit=gate_circuit, inv=False)

    _build_P_gates(gate=gate, nqubits=nqubits, circuit=gate_circuit, inv=True)

    _build_Q_gates(nqubits=nqubits-1, circuit=gate_circuit, inv=True)

    circuit.compose(gate_circuit, controls + [target], inplace=True)