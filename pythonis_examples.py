
#######
#
# Testing scripts
#
#######


# test line to use on my own setup - gene network generated from TDCor inference - around 250 genes
a, b, c = ImportBooleanModel("C:\Work\Datasets\LRPNetwork-parallel-20150525", "l_gnp-190515.txt",
                             "TDCor6.32_output_Missinglink_parallel-190515-cytoscape.txt")
flow, stable, start = RunBooleanModel(a, b, 't', 't', initial_state_number=5, initial_state_choice='random',
                                      stimulus='constant')

# Script to run the algo a bunch of time and plot the resulting scatter basin size = f(stable state size)
IterPlot(a, b, 50, 5)

# Script to check the activity of a given gene in all stable states
CheckGeneStableStates('SHR', a, stable)

# Filter through all the stable states to get the lists of genes which change activity
FilterStableStates(a, stable)

# Simple gene networks for testing purpose :
# one 2 states cycle, basin 2 ; one 6 states cycle, basin 6
a1 = ['A', 'B', 'C']
b1 = {'A': [['1', 'B']], 'B': [['1', 'C']], 'C': [['-1', 'A']]}
flow1, stable1, start1 = RunBooleanModel(a1, b1, 't', 't', 'all')

# one 3 states cycle, basin 6 ; one stable state, basin 10
a2 = ['A', 'B', 'C', 'D']
b2 = {'A': [['-1', 'B']], 'B': [['-1', 'C']], 'C': [['1', 'A'], ['1', 'B'], ['-1', 'D']], 'D': [['1', 'A']]}
flow2, stable2, start2 = RunBooleanModel(a2, b2, 't', 't', 'all')

# 2 stable state, basin 4 and 12
a3 = ['A', 'B', 'C', 'D']
b3 = {'A': [['-1', 'D']], 'B': [['-1', 'C']], 'C': [['1', 'A'], ['1', 'C'], ['-1', 'D']], 'D': [['1', 'A']]}
flow3, stable3, start3 = RunBooleanModel(a3, b3, 't', 't', 'all')

# Simple tests
a4 = ['A', 'B']
b4 = {'A': [['1', 'B']]}
flow4, stable4, start4 = RunBooleanModel(a4, b4, 't', 't', 'all')

a5 = ['A', 'B']
b5 = {'A': [['-1', 'B']]}
flow5, stable5, start5 = RunBooleanModel(a5, b5, 't', 't', 'all')

a6 = ['A', 'B']
b6 = {'A': [['1', 'B'], ['1', 'A']]}
flow6, stable6, start6 = RunBooleanModel(a6, b6, 't', 't', 'all')

a7 = ['A', 'B']
b7 = {'A': [['-1', 'B']], 'B': [['1', 'B']]}
flow7, stable7, start7 = RunBooleanModel(a7, b7, 't', 't', 'all')

a8 = ['A', 'B']
b8 = {'A': [['1', 'B'], ['-1', 'A']]}
flow8, stable8, start8 = RunBooleanModel(a8, b8, 't', 't', 'all')

a9 = ['A']
b9 = {'A': [['1', 'A']]}
flow9, stable9, start9 = RunBooleanModel(a9, b9, 't', 't', 'all')

a10 = ['A']
b10 = {'A': [['-1', 'A']]}
flow10, stable10, start10 = RunBooleanModel(a10, b10, 't', 't', 'all')

# Complex tests - just changing the nature of the interactions, not the topology - checking if we can generate different combinations of stable states and basin size

# 1 stable cycle (size 2), basin size 128
a11 = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
b11 = {'A': [['-1', 'B'], ['-1', 'G']], 'B': [['1', 'C'], ['1', 'A'], ['1', 'F']],
       'C': [['1', 'A'], ['-1', 'C'], ['1', 'D']], 'D': [['-1', 'A'], ['1', 'B']], 'E': [['1', 'C'], ['-1', 'A']],
       'F': [['1', 'C'], ['-1', 'E']], 'G': [['1', 'D']]}
flow11, stable11, start11 = RunBooleanModel(a11, b11, 't', 't', 'all')

# 3 stable states (1, 8) - (4, 96) - (2, 24)
a12 = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
b12 = {'A': [['1', 'B'], ['-1', 'G']], 'B': [['1', 'C'], ['-1', 'A'], ['1', 'F']],
       'C': [['1', 'A'], ['1', 'C'], ['1', 'D']], 'D': [['1', 'A'], ['1', 'B']], 'E': [['1', 'C'], ['1', 'A']],
       'F': [['1', 'C'], ['1', 'E']], 'G': [['1', 'D']]}
flow12, stable12, start12 = RunBooleanModel(a12, b12, 't', 't', 'all')

# 3 stable states (3, 96) - (2, 24) - (2, 8)
a13 = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
b13 = {'A': [['1', 'B'], ['1', 'G']], 'B': [['1', 'C'], ['-1', 'A'], ['-1', 'F']],
       'C': [['1', 'A'], ['-1', 'C'], ['1', 'D']], 'D': [['1', 'A'], ['-1', 'B']], 'E': [['1', 'C'], ['1', 'A']],
       'F': [['1', 'C'], ['1', 'E']], 'G': [['1', 'D']]}
flow13, stable13, start13 = RunBooleanModel(a13, b13, 't', 't', 'all')

# for network pruning test 
pruning_test = {'Z1':[['1','Z2']],'Z2':[['1','Z3']],'Z3':[['1','A']], 'A': [['-1', 'B'], ['-1', 'C']], 'B': [['1', 'D'],['1','J']],
           'C': [['1', 'E']], 'D': [['-1', 'B'], ['1', 'F']], 'E': [['1', 'C'], ['-1', 'G']],
           'F': [['1', 'H']], 'G': [['1', 'I']]}
core_test = ExtractCoreNetwork(pruning_test)
print pruning_test
print core_test

