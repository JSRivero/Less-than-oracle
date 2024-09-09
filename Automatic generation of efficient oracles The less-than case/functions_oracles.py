# Imports
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import Diagonal
from linear_mc_gate import linear_mc_gate

import numpy as np


#######################
###### UTILITIES ######
#######################

def to_binary(number, nbits=None):
    
    '''
    This fucntion transforms an integer to its binary form (string).
    If a determined number of bits is required (more than the needed ones),
    it can be passed as a parameter too, nbits, None by default.
    It is needed that the number of bits passed as a parameter is larger
    than the number of bits needed to write the number in binary. 

    Input:
    number: integer (int).
    nbits: integer (int), None by default

    Output:
    binary: string (str) containing the number in its binary form.
    It writes 0s in front if nbits is larger than the number of bits needed
    to write the binary form.
    '''

    if nbits is None:
        return bin(number)[2:]
    else:
        binary = bin(number)[2:]
        if nbits < len(binary):
            print('Error, nbits must be larger than %d.'%(len(binary)))
        else:
            return '0' * (nbits - len(binary)) + binary


def multi_control_z(nqubits):
    '''
    Function to create a multi-controlled Z gate.

    Input:
    nqubits: Integer (int) of the number of qubits in the gate (controls and target)
       This means that the gate has nqubits-1 controls and 1 target.

    Output:
    circuit: QuantumCircuit containing a multi-controlled Z gate.
      It has to be transformed with method .to_gate() to append to a QuantumCircuit larger.

    Example:

    main_circuit = QuantumCircuit(nqubits)

    gate_multi_z = multi_control_z(nqubits)

    main_circuit.append(gate_multi_z.to_gate(), range(nqubits))
    '''
    circuit=QuantumCircuit(nqubits,name=' CZ (%d)' %(nqubits))
    # circuit.h(nqubits-1)
    # gate = MCXGate(nqubits-1)
    # circuit.append(gate, range(nqubits))
    # circuit.h(nqubits-1)
    # return circuit

    xgate = np.array([[1,0],[0,-1]])
    controls = list(range(nqubits-1))
    target = nqubits-1
    linear_mc_gate(gate=xgate, circuit=circuit, controls=controls, targ=target)
    return circuit


#################################################
############ ORACLE LESS-THAN ORACLE ############
#################################################


def oracle_less_than(number, nqubits, name=None):

    '''
    This function builds a quantum circuit, an oracle, which marks with a pi-phase
    those states which represent numbers strictly smaler than the number given by parameter.

    The procedure is almost the same for all numbers, with the only exception of a difference
    if the first bit of the number in binary is 1 or 0.

    Input:
    number: integer (int) containing the objective number,
       or a string (str) with the binary representation of such number.
    nqubits: integer (int) number of qubits of the circuit.
       It must be larger than the number of digits of the binary representation of number.
    name: string (str), default None, name of the circuit.

    Output:
    circuit: QuantumCircuit which marks with phase pi the states which
    represent in binary the numbers strictly smaller than number.
    '''

    # Construction of the circuit
    if name:# If name is provided give such name to the circuit
        circuit = QuantumCircuit(nqubits, name=name)
    else: # Otherwise, the name is just " < number"
        circuit = QuantumCircuit(nqubits, name = ' < %d '%number)

    # Binary representation of the number
    num_binary = to_binary(number, nqubits)
    
    # Discard the 0s at the end, as they will not be used and save
    # unnecessary X gates
    num_binary = num_binary.rstrip('0')

    
    if num_binary[0] == '1':
        # If the first digit is 1
        # Mark all the states of the form |0q1...>
        circuit.x(nqubits-1)
        circuit.z(nqubits-1)
        circuit.x(nqubits-1)
    else:
        # If first digit is 0
        # Apply X gate to first qubit
        circuit.x(nqubits-1)
    
    # For loop on the remaining digits
    for position1, value in enumerate(num_binary[1:]):
        # Rename the position as it starts with 0 in the second bit and
        # we want it to be 1.
        position = position1 + 1

        if value == '0':
            # If the digit is 0
            # Just apply a X gate
            circuit.x(nqubits-position-1)
        else:
            # If the digit bi is 1
            # Apply a multi-controlled Z gate to mark states of the shape:
            # |bn...bi+1 0 qi-1...q1>
            # where bn,...,bi+1 are the first n-i bits of m, which is of the shape bn...bi+1 1 bi-1...b1
            # because we just checked that bi is 1.
            # Hence, the numbers of the form bn...bi+1 0 qi-1...q1 are smaller than m.
            circuit.x(nqubits-position-1)
            multi_z = multi_control_z(position + 1)
            circuit.append(multi_z.to_gate(), range(nqubits-1, nqubits-position-2, -1))
            circuit.x(nqubits-position-1)
    
    for position, value in enumerate(num_binary):
        # Apply X gates to qubits in position of bits with a 0 value
        if value == '0':
            circuit.x(nqubits-position-1)
        else:
            pass
    
    return circuit


#################################################
###### ORACLES BUILT WITH LESS-THAN ORACLE ######
#################################################


def oracle_greater_than(number, nqubits, name=None):

    '''
    This function builds a quantum circuit, an oracle, which marks with a pi-phase
    those states which represent numbers strictly larger than the number given by parameter.

    The procedure uses the function oracle_less_than and introduces a global phase.

    Input:
    number: integer (int) containing the objective number,
       or a string (str) with the binary representation of such number.
    nqubits: integer (int) number of qubits of the circuit.
       It must be larger than the number of digits of the binary representation of number.
    name: string (str), default None, name of the circuit.

    Output:
    circuit: QuantumCircuit which marks with phase pi the states which
    represent in binary the numbers strictly greater than number.
    '''
    if name is None:
        name = ' > %d'%(number)
    else:
        name = name
    
    circuit = oracle_less_than(number=number+1, nqubits=nqubits, name=name)

    circuit.z(0)
    circuit.x(0)
    circuit.z(0)
    circuit.x(0)
    
    return circuit


def oracle_interval(lower_boundary, upper_boundary, nqubits, name=None):

    '''
    This function builds a quantum circuit, an oracle, which marks with a pi-phase
    those states which represent numbers in the interval (lower_boundary, upper_boundary).

    The procedure uses the function oracle_less_than and oracle_greater_than.

    Input:
    lower_boundary: integer (int) containing the lower boundary,
       or a string (str) with the binary representation of such number.
    upper_boundary: integer (int) containing the upper boundary,
       or a string (str) with the binary representation of such number.
    nqubits: integer (int) number of qubits of the circuit.
       It must be larger than the number of digits of the binary representation of number.
    name: string (str), default None, name of the circuit.

    Output:
    circuit: QuantumCircuit which marks with phase pi the states which
    represent in binary the numbers in the interval (lower_boundary, upper_boundary).
    '''

    if name is None:
        name = '(%d,%d)'%(lower_boundary, upper_boundary)
    else:
        name = name

    circuit = QuantumCircuit(nqubits, name=name)

    circuit.append(oracle_less_than(number=upper_boundary, nqubits=nqubits), range(nqubits))
    circuit.append(oracle_greater_than(number=lower_boundary, nqubits=nqubits), range(nqubits))
    
    circuit.z(0)
    circuit.x(0)
    circuit.z(0)
    circuit.x(0)

    return circuit


#############################################
###### ORACLES BUILT WITH DIAONAL GATE ######
#############################################

def array_less_than(number:int, N:int):
    '''
    Creates a 1-D vector of N elements where the first "number" elements
    are -1s and the rest are 1s.
    '''
    v = np.array([-1]*number + [1]*(N-number))
    return v

def diagonal_oracle_less_than(number:int, nqubits:int, name=None):
    '''
    This function builds a quantum circuit, an oracle, which marks with a pi-phase
    those states which represent numbers strictly smaler than the number given by parameter
    by using the method Diagonal from Qiskit.

    Input:
    number: integer (int) containing the objective number,
       or a string (str) with the binary representation of such number.
    nqubits: integer (int) number of qubits of the circuit.
       It must be larger than the number of digits of the binary representation of number.
    name: string (str), default None, name of the circuit.

    Output:
    circuit: QuantumCircuit which marks with phase pi the states which
    represent in binary the numbers strictly smaller than number.
    '''

    v = array_less_than(number=number, N=2**nqubits)

    # Build diagonal circuit
    diag_gate = Diagonal(v)
    diagonal = QuantumCircuit(nqubits, name=name)
    diagonal.append(diag_gate, range(nqubits))

    return diagonal



################################################
###### CALCULATE DEPTH WITH GIVEN BACKEND ######
################################################

def decompose_circuit(circuit:QuantumCircuit, backend, reps:int=50, opt_level=3):
    '''
    Decomposes a circuit according to a backend.

    Input:
        - circuit (QuantumCircuit): circuit to be decomposed

        - backend: backend to use. It may be a str (when statevector), or a backend object

        - reps (int): number of repetitions in decompose. By default uses 50 which is
        usually enough to decompose any circuit.

        - opt_level (int in 0, 1, 2, 3): optimization level of the circuit. By default 3 to
        obtain the maximum optimization.
    '''

    circuit_transpiled = transpile(circuits=circuit, backend=backend, optimization_level=opt_level)
    circuit_decomposed = circuit_transpiled.decompose(reps=reps)
    
    return circuit_decomposed