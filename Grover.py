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
Apply Hadamards to all specified qubits (if apply_to_list[index] == 1)
Designed for a large amount of Hadamards being applied at once
"""
def apply_H(program, apply_to_list):
    for index,qubit in enumerate(apply_to_list):
        if qubit == 1:
            program += H(index)
            
            
    
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
    
        
def get_Z0(n):
#Apply the negative in G here
    gate = np.zeros((2**n, 2**n), dtype=int)
    gate[0][0] = 1
    for i in range(1, 2**n):
        gate[i][i] = -1
#    for i in range(2**n):
#        print(gate[i])
    return DefGate("Z0", gate)
    

def get_Zf(f, n):
    bitstrings = list()
    get_bitstring_permutations(0, bitstrings, n, [0]*n)
    gate = np.zeros((2**n, 2**n), dtype=int)
    for i in range(2**n):
        gate[i][i] = -1 if f(bitstrings[i]) == 1 else 1
    return DefGate("Zf", gate)
    
#Assume n <= 9
def grovers_algorithm(f, n):
    #initialize circuit
    program = initialize([0]*n)
    qubits = list(range(n))
    z0_gate = get_Z0(n)
    zf_gate = get_Zf(f,n)
    Z0 = z0_gate.get_constructor()
    Zf = zf_gate.get_constructor()
    iteration_count = floor(pi/4 * sqrt(n))
    h_qubits = [1]*n
    for i in range(iteration_count):
        #Apply Zf
        program += zf_gate
        program += Zf(*qubits)
        #Apply H to all qubits
        apply_H(program, h_qubits)
        #Apply -Z0
        program += z0_gate
        program += Z0(*qubits)
        #Apply H to all qubits
        apply_H(program, h_qubits)
    #Measure
    #print(program)
    with local_forest_runtime():
       qvm = get_qc('9q-square-qvm')
       results = qvm.run_and_measure(program, trials=10)
       #print(results[0])
       for i in range(len(results)):
#           for item in results:
#                print(type(item))
           new = list()
           for j in range(len(results[i])):
              new.append(results[i][j])
           if f(new) == 1: return 1
       return 0
        
    
def f(list):
    for l in list:
        if l == 1:
            return 0
    return 0

print(grovers_algorithm(f, 5))


