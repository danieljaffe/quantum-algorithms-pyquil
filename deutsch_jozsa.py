from pyquil import Program, get_qc
from pyquil.gates import *
from pyquil.quilbase import DefPermutationGate
from pyquil.quil import DefGate
from pyquil.api import local_forest_runtime
from math import floor, pi, sqrt
import numpy as np


def get_bitstring_permutations(index, lst, n, args):
    """
    This function populates a list with all the combinations of bit strings of length n
    """
    if (index == n):
        # append combination to list
        lst.append(list.copy(args))
    else:
        # handle the case where the preceding bit is 0
        args[index] = 0
        get_bitstring_permutations(index + 1, lst, n, args)

        # handle the case where the preceding bit is 1
        args[index] = 1
        get_bitstring_permutations(index + 1, lst, n, args)


def generate_uf(f, n, name):
    """
    Parameters: f is an anonymous function and n is the number of bits in input: f:{0,1}^n -> {0,1}
    This function returns an oracle gate representing the function f
        for all x in {0,1}^n and y in {0,1}, the desired result is of the oracle is mapping the input
        |x,y> to |x, y + f(x)> where + is addition modulo 2. The function first finds the list of bitstring
        permutations of n bits, it then establishes a mapping which is representative of the decimal number of the
        bitstring represents. It then determines for each |x,y>, it calculates |x, f(x) + y>. Finally it constructs a
        permutation gate which treats each permutation as a different basis vector in the 2^(n+1) dimensional complex
        hilbert space that represents a system of n + 1 qubits. The permutation gate is returned.
    """
    # generate list of all bitstrings of size n
    bitstrings = list()
    get_bitstring_permutations(0, bitstrings, n + 1, [0] * (n + 1))
    # initialize mapping and permutation list
    perm_dict = dict()
    perm_list = list()
    # populate mapping
    for permutation, bitstring in enumerate(bitstrings):
        perm_dict["".join(str(bit) for bit in bitstring)] = permutation
    # Send each |xy> to |x, f(x) + y>
    for bitstring in bitstrings:
        params = bitstring[:n]
        params.append((f(params) + bitstring[-1]) % 2)
        perm_list.append(perm_dict["".join(str(bit) for bit in params)])
    # Create and return permutation gate
    return DefPermutationGate(name, perm_list)


def apply_H(program, apply_to_list):
    """
    Apply Hadamards to all specified qubits (if apply_to_list[index] == 1).
    Designed for a large amount of Hadamards being applied at once.
    """
    for index, qubit in enumerate(apply_to_list):
        if qubit == 1:
            program += H(index)
    return program


def initialize(states):
    """
    This function sets initial states and applies Hadamards to each qubit
    Note: apply_H isn't called because it is actually more efficient to initialize in one loop as opposed to 2.
    """
    program = Program()
    for index, state in enumerate(states):
        if state == 1:
            program += X(index)
        program += H(index)
    return program


def deutsch_jozsa_algorithm(f, n):
    """
    This function is intended to determine if f is constant or balanced for a given function f s.t.
        f:{0,1}^n -> {0,1}. The algorithm initializes the qubits with H for the first n qubits and X and H for the last
        qubit. The algorithm then constructs a Uf oracle gate based on the function input f. It then applies Uf to all
        the qubits and applies H to the first n qubits. Finally, the simulator is run on the circuit and measures the
        results. If upon measurement, the first n qubits are all 0, 1 is returned and the function is constant,
        otherwise 0 is returned and the function is balanced.
    This function has an anonymous function and integer n as parameters.
    This function uses 9q-squared-qvm, so it assumes that n <= 9.
    """
    # apply H to first n qubits and X H to the last qubit (ancilla qubit)
    initialize_list = [0] * n
    initialize_list.append(1)
    qubits = list(range(len(initialize_list)))
    program = initialize(initialize_list)
    # Generate Uf oracle from f (anonymous function)
    uf_gate = generate_uf(f, n, "Uf")
    Uf = uf_gate.get_constructor()
    # Add Uf to circuit
    program += uf_gate
    program += Uf(*qubits)
    # Apply Hadamards to first n qubits
    apply_to_list = [1] * n
    apply_to_list.append(0)
    program = apply_H(program, apply_to_list)
    # print(program)
    # Run simulator
    with local_forest_runtime():
        # Assume n <= 9
        qvm = get_qc('9q-square-qvm')
        # 1 trial because DJ is deterministic
        results = qvm.run_and_measure(program, trials=1)
        # check if first n qubits are all 0
        for i in range(n):
            if (results[i] != 0):
                return 0
        return 1
    # in case of failure, return 0
    return 0


def f(args):
    return 1
    # return (args[0] + args[1])%2


print(deutsch_josza_algorithm(f, 2))
