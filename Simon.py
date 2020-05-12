from pyquil import Program, get_qc
from pyquil.gates import *
from pyquil.quilbase import DefPermutationGate
from pyquil.quil import DefGate
from pyquil.api import local_forest_runtime
from math import *
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
    This function sets initial states and applies Hadamards to each qubit.
    Note: apply_H isn't called because it is actually more efficient to initialize in one loop as opposed to 2.
    """
    program = Program()
    for index, state in enumerate(states):
        if state == 1:
            program += X(index)
        program += H(index)
    return program


def generate_uf_simons(f, n, name):
    """
    Parameters: f is an anonymous function and n is the number of bits in input: f:{0,1}^n -> {0,1}^n
    This function returns an oracle gate representing the function f for all x in {0,1}^n and y in {0,1}^n,
    the desired result is of the oracle is mapping the input x,y> to |x, y + f(x)> where + is addition modulo 2.
    The function first finds the list of bitstring permutations of n bits, it then establishes a mapping which is
    representative of the decimal number of the bitstring represents. For each |x,y>, it calculates
    |x, y + f(x)>. Finally it constructs a permutation gate which treats each permutation as a different basis vector
    in the 2^(n+1) dimensional complex Hilbert space that represents a system of 2*n qubits.
    Returns: the permutation gate
    """

    # generate list of all bitstrings of size n
    bitstrings = list()
    get_bitstring_permutations(0, bitstrings, 2 * n, [0] * (2 * n))

    # initialize mapping and permutation list
    perm_dict = dict()
    perm_list = list()

    # populate mapping
    for permutation, bitstring in enumerate(bitstrings):
        perm_dict["".join(str(bit) for bit in bitstring)] = permutation

    # Send each |xy> to |x, f(x) + y>
    for bitstring in bitstrings:
        params = bitstring[:n]
        params2 = bitstring[n:2 * n]
        f_values = f(bitstring[:n])
        for i in range(n):
            params.append((params2[i] + f_values[i]) % 2)
        perm_list.append(perm_dict["".join(str(bit) for bit in params)])

    # Create and return permutation gate
    return DefPermutationGate(name, perm_list)


def simons_solver(Y, n):
    """
    Inputs: Y is a linear system of n-1 equations in matrix form. n is the dimension of the input into f.
    This function acts as a binary linear matrix solver.
    Returns: the key string s, if found, or the zero bitstring
    """
    # Create all possible bit strings to test for s
    bitstrings = []
    get_bitstring_permutations(0, bitstrings, n, [0] * n)

    # For each possible s, test to see if it's a candidate
    for s in bitstrings:
        if s == [0] * n:
            continue
        candidate = True

        # For each equation in Y, bit by bit test that y*s = [0]*n
        for y in Y:
            value = 0
            for i in np.arange(n):
                value = value + s[i] * y[i]

            # If a bit doesn't evaluate to 0...
            if (value % 2 == 1):
                candidate = False
        if (candidate):
            return s

    return [0] * n


def simons_algorithm(f, n):
    """
    Inputs: f is a blackbox function (f:{0,1}^n -> {0,1}^n) that is either one-to-one or two-to-onw. n is the
    dimension of the input into f. This function finds and returns the key s, if one exists, for a two=to-one
    function by first creating a matrix U_f that represents f, then applying the appropriate quantum gates to
    generate a linear equation. By running the circuit n-1 times, we generate n-1 equations that we then feed into a
    classical solver. The Classical solver returns s.
    Returns: the key string s, if found, or the zero bitstring
    """
    # Initialize the program and apply the Hadamard gate to the first n qubits.
    program = initialize([0] * n)

    # Initialize an extra n bits by applying the Identity gate.
    for index in range(n, (2 * n)):
        program += I(index)
    # Generate the U_f gate for f
    uf_simons_gate = generate_uf_simons(f, n, "Uf_Simons")
    Uf_simons = uf_simons_gate.get_constructor()

    # Initialize qubits and apply the U_f gate to all qubits
    qubits = list(range(2 * n))
    program += uf_simons_gate
    program += Uf_simons(*qubits)

    # Pick the qubits to which to apply the final Hadamard gates.
    apply_to_list = [1] * n
    np.append(apply_to_list, [0] * n)
    program = apply_H(program, apply_to_list)

    with local_forest_runtime():
        qvm = get_qc('9q-square-qvm')
        s_trials = list()

        # Conduct 20 trials of the algorithm - it's not deterministic
        for i in range(20):

            # Run the circuit n-1 times to generate n-1 equations
            measurements = qvm.run_and_measure(program, trials=n - 1)
            Y = np.zeros((n - 1, n))

            # Parse the resulting relevant measurements from the circuit and store in a matrix
            for index in range(n):
                for j in range(len(measurements[index])):
                    Y[j][index] = measurements[index][j]

            # Plug Y into the binary matrix solver and store in a list of possible s strings (s_trials)
            s_primes = simons_solver(Y, n)
            s_trials.append(s_primes)

        # Check if each s satisfies the condition needed to be a 'key' and return it if the condition is met
        for trial in s_trials:
            if (f([0] * n) == f(trial)):
                return trial

        # If no such s is found, return the 0 bitstring to indicate s doesn not exist and the function is one-to-one
        return [0] * n

# def f(x):
#     if x == [0, 0, 0]:
#         return [1, 0, 1]
#     elif x == [0, 0, 1]:
#         return [0, 1, 0]
#     elif x == [0, 1, 0]:
#         return [0, 0, 0]
#     elif x == [0, 1, 1]:
#         return [1, 1, 0]
#     elif x == [1, 0, 0]:
#         return [0, 0, 0]
#     elif x == [1, 0, 1]:
#         return [1, 1, 0]
#     elif x == [1, 1, 0]:
#         return [1, 0, 1]
#     elif x == [1, 1, 1]:
#         return [0, 1, 0]
#     else:
#         return NULL
#
#
# def f2(x):
#     if x == [0, 0]:
#         return [0, 1]
#     elif x == [0, 1]:
#         return [1, 1]
#     elif x == [1, 0]:
#         return [0, 1]
#     elif x == [1, 1]:
#         return [1, 1]
#
#
# def f3(x):
#     if x == [0]:
#         return [1]
#     elif x == [1]:
#         return [1]
#
#
# s = simons_algorithm(f2, 2)
# print(s)
