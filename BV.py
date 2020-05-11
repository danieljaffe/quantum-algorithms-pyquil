from pyquil import Program, get_qc
from pyquil.gates import *
from pyquil.quilbase import DefPermutationGate
from pyquil.quil import DefGate
from pyquil.api import local_forest_runtime
from math import floor, pi, sqrt
import numpy as np

#This function populates a list with all the combinations of bit strings of length n
def get_bitstring_permutations(index, lst, n, args):
    if(index == n):
        #append combination to list
        lst.append(list.copy(args))
    else:
        #handle the case where the preceding bit is 0
        args[index] = 0
        get_bitstring_permutations(index + 1, lst, n, args)
        
        #handle the case where the preceding bit is 1
        args[index] = 1
        get_bitstring_permutations(index + 1, lst, n, args)

"""
Parameters: f is an anoynmous function and n is the number of bits in input: f:{0,1}^n -> {0,1}
This function returns an oracle gate representing the function f
    for all x in {0,1}^n and y in {0,1}, the desired result is of the oracle is mapping the input |x,y> to |x, y + f(x)> where + is addition modulo 2
    The function first finds the list of bitstring permutations of n bits, it then establishes a mapping which is representative of the decimal number of the bitstring reprents. It then determines for each |x,y>, it calculates |x, f(x) + y>. Finally it constructs a perumtation gate which treats each permutation as a different basis vector in the 2^(n+1) dimensional complex hilbert space that represents a system of n + 1 qubits
    The permutation gate is returned.
"""
def generate_uf(f,n,name):
    #generate list of all bitstrings of size n
    bitstrings = list()
    get_bitstring_permutations(0, bitstrings, n+1, [0]*(n+1))
    #initialize mapping and permutation list
    perm_dict = dict()
    perm_list = list()
    #populate mapping
    for permutation, bitstring in enumerate(bitstrings):
        perm_dict["".join(str(bit) for bit in bitstring)] = permutation
    #Send each |xy> to |x, f(x) + y>
    for bitstring in bitstrings:
        params = bitstring[:n]
        params.append((f(params) + bitstring[-1])%2)
        perm_list.append(perm_dict["".join(str(bit) for bit in params)])

    #Create and return permutation gate
    return DefPermutationGate(name, perm_list)
    
"""
Apply Hadamards to all specified qubits (if apply_to_list[index] == 1)
Designed for a large amount of Hadamards being applied at once
"""
def apply_H(program, apply_to_list):
    for index,qubit in enumerate(apply_to_list):
        if qubit == 1:
            program += H(index)
    return program
        

"""
This function sets initial states and applies Hadamards to each qubit
Note: apply_H isn't called because it is actually more efficient to initialize in one loop as opposed to 2.
"""
def initialize(states):
    program = Program()
    for index, state in enumerate(states):
        if state == 1:
            program += X(index)
        program += H(index)
    return program
    
    
    
def bernstein_vazirani_algorithm(f,n):
    initialize_list = [0]*n
    b = f(initialize_list)
    initialize_list.append(1)
    qubits = list(range(len(initialize_list)))
    program = initialize(initialize_list)
    uf_gate = generate_uf(f,n,"Uf")
    Uf = uf_gate.get_constructor()
    program += uf_gate
    program += Uf(*qubits)
    apply_to_list = [1]*n
    apply_to_list.append(0)
    program = apply_H(program, apply_to_list)
    print(program)
    with local_forest_runtime():
        qvm = get_qc('9q-square-qvm')
        #1 trial because DJ is deterministic
        results = qvm.run_and_measure(program, trials=1)
        ones = 0
       # print(results)
        a = []
        for i in range(n):
            a.append(results[i][0])
        return a,b
    #in case of failure
    return b, None
    #return (args[0] + args[1])%2
    
def f(args):
    return args[1]

a,b = bernstein_vazirani_algorithm(f, 2)
print(a, b)

