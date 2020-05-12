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


def apply_H(program, apply_to_list):
    """
    Apply Hadamards to all specified qubits (if apply_to_list[index] == 1)
    Designed for a large amount of Hadamards being applied at once
    """
    for index, qubit in enumerate(apply_to_list):
        if qubit == 1:
            program += H(index)


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


def get_Z0(n):
    """
    This function generates the Z0 gate satisfying the conditions for x in {0,1}^n Z0|x> -> -|x> iff x = 0^n
        otherwise Z0|x> -> |x>
    The parameter to this function is only the size n, a 2^n x 2^n dimensional matrix is created satisfying the
        conditions above.
    This function has one dependency, the DefGate function defined in pyquil.quil
    This function is designed to absorb the negative in G, so the returned gate is actually -Z0
    Returns -Z0
    """
    # Create a 2^n x 2^n matrix with all 0's
    gate = np.zeros((2 ** n, 2 ** n), dtype=int)
    # since this is -Z0, set first element to 1 not -1
    gate[0][0] = 1
    # set all other elements on the diagonal to -1, again not 1 because this is -Z0
    for i in range(1, 2 ** n):
        gate[i][i] = -1
    # Return gate
    return DefGate("Z0", gate)


def get_Zf(f, n):
    """
    This function generates the Zf gate satisfying the condition for x in {0,1}^n where Zf|x> -> (-1)^f(X)|x>
    This function requires that f(x) be calculated for all x, so f is passed as an anonymous function, the other
        parameter is n.
    The function has one dependency, the DefGate function defined in pyquil.quil
    This function finds all permutations of bitstrings of length n, then initializes a 2^n x 2^n matrix of all 0's,
        and sets all elements along the diagonal to either 1 or -1 depending on f(x)
    Finally a gate representation of this matrix is returned.
    """
    # generate bitstring permutations
    bitstrings = list()
    get_bitstring_permutations(0, bitstrings, n, [0] * n)
    # initialize a 2^n x 2^n matrix of all 0's
    gate = np.zeros((2 ** n, 2 ** n), dtype=int)
    # set diagonals of matrix based on f(x)
    for i in range(2 ** n):
        gate[i][i] = -1 if f(bitstrings[i]) == 1 else 1
    # create and return gate
    return DefGate("Zf", gate)


def grovers_algorithm(f, n):
    """
    This function is intended to determine if there exists an x in {0,1}^n s.t. f(x) = 1 for a given function f s.t.
        f:{0,1}^n -> {0,1}. The algorithm first constructs Zf, -Z0 gates, initializes with Hanamard matrices, and
        applies G = -H^n o Z0 o H^n o Zf. This algorithm is not deterministic, so G is applied multiple times. More
        specifically, G is run (pi / 4 * sqrt(n)) times. Furthermore, there are 10 trials to minimize the chance of a
        false negative.
    This function has an anonymous function and integer n as parameters.
    This function runs the algorithm as described for each 10 trials, and then checks if for any of the outputted states
        x, if f(x) = 1. If this is true, then 1 is returned, otherwise 0 is returned. The function returns 0 if there
        is an issue with the simulator.
    This function uses 9q-squared-qvm, so it assumes that n <= 9
    """
    # Apply Hadamards to all qubits
    program = initialize([0] * n)
    qubits = list(range(n))
    # Define and generate Z0 gate (really -Z0)
    z0_gate = get_Z0(n)
    Z0 = z0_gate.get_constructor()
    # Define and generate Zf gate
    zf_gate = get_Zf(f, n)
    Zf = zf_gate.get_constructor()
    # Determine the number of times to apply G
    iteration_count = floor(pi / 4 * sqrt(n))
    h_qubits = [1] * n
    # Apply G iteration_count times
    for i in range(iteration_count):
        # Apply Zf
        program += zf_gate
        program += Zf(*qubits)
        # Apply H to all qubits
        apply_H(program, h_qubits)
        # Apply -Z0
        program += z0_gate
        program += Z0(*qubits)
        # Apply H to all qubits
        apply_H(program, h_qubits)
    # Measure
    # Try to run simulator
    with local_forest_runtime():
        # assumes that n <= 9
        qvm = get_qc('9q-square-qvm')
        # run circuit and measure the qubits, 10 trials
        results = qvm.run_and_measure(program, trials=10)
        # Iterate through different results of the different trials and check if f(x) = 1
        for i in range(len(results)):
            new = list()
            for j in range(len(results[i])):
                new.append(results[i][j])
            if f(new) == 1: return 1
        # return 0 if simulator fails
        return 0
