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
    This function sets initial states and applies Hadamards to each qubit.
    Note: apply_H isn't called because it is actually more efficient to initialize in one loop as opposed to 2.
    """
    program = Program()
    for index, state in enumerate(states):
        if state == 1:
            program += X(index)
        program += H(index)
    return program

def generate_uf_simons(f,n,name):
    #generate list of all bitstrings of size n
    bitstrings = list()
    get_bitstring_permutations(0, bitstrings, 2*n, [0]*(2*n))
    #initialize mapping and permutation list
    perm_dict = dict()
    perm_list = list()
    #populate mapping
    for permutation, bitstring in enumerate(bitstrings):
        perm_dict["".join(str(bit) for bit in bitstring)] = permutation
    #Send each |xy> to |x, f(x) + y>
    for bitstring in bitstrings:
        params = bitstring[:n]
        params2 = bitstring[n:2*n]
        f_values = f(bitstring[:n])
        for i in range(n):
            params.append((params2[i]+f_values[i])%2)
        perm_list.append(perm_dict["".join(str(bit) for bit in params)])

    #Create and return permutation gate
    return DefPermutationGate(name, perm_list)

def simons_solver(Y, n):
    bitstrings = []
    get_bitstring_permutations(0, bitstrings, n, [0]*n)
    for s in bitstrings:
        if s == [0]*n:
            continue
        candidate = True
        for y in Y:
            value = 0
            for i in np.arange(n):
                value = value + s[i]*y[i]
            if(value%2 == 1):
                candidate = False
        if (candidate):
            return s

    return [0]*n

def simons_algorithm(f, n):
    program = initialize([0]*n)
    for index in range(n, (2*n)):
        program += I(index)
    uf_simons_gate = generate_uf_simons(f,n,"Uf_Simons")
    Uf_simons = uf_simons_gate.get_constructor()
    qubits = list(range(2*n))
    program += uf_simons_gate
    program += Uf_simons(*qubits)
    apply_to_list = [1]*n
    np.append(apply_to_list, [0]*n)
    program = apply_H(program, apply_to_list)

    with local_forest_runtime():
        qvm = get_qc('9q-square-qvm')
        s_trials = list()
        for i in range(20):
            measurements = qvm.run_and_measure(program, trials=n-1)
            Y = np.zeros((n-1, n))
            for index in range(n):
                for j in range(len(measurements[index])):
                    Y[j][index] = measurements[index][j]
            s_primes = simons_solver(Y, n)
            s_trials.append(s_primes)
        for trial in s_trials:
            if (f([0]*n) == f(trial)):
                return trial
        return [0]*n

def f(x):
    if x == [0,0,0]:
        return [1,0,1]
    elif x == [0,0,1]:
        return [0,1,0]
    elif x == [0,1,0]:
        return [0,0,0]
    elif x == [0,1,1]:
        return [1,1,0]
    elif x == [1,0,0]:
        return [0,0,0]
    elif x == [1,0,1]:
        return [1,1,0]
    elif x == [1,1,0]:
        return [1,0,1]
    elif x == [1,1,1]:
        return [0,1,0]
    else:
        return NULL

def f2(x):
    if x == [0,0]:
        return [0,1]
    elif x == [0,1]:
        return [1,1]
    elif x == [1,0]:
        return [0,1]
    elif x == [1,1]:
        return [1,1]

def f3(x):
    if x == [0]:
        return [1]
    elif x == [1]:
        return [1]

s = simons_algorithm(f2, 2)
print(s)
