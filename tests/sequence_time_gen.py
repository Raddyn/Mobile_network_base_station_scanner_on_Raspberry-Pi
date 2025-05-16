import numpy as np
import time
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import core.sequence_generators as seq_gen

'''
Testing of the sequence generation time
for PSS and SSS.
'''

# PSS
#------------------------------------------------------------------

start_time = time.time()
pss = np.zeros((3, 62), dtype=complex)
for i in range(3):
    pss[i,:] = seq_gen.pss_gen(i)

# SSS
#------------------------------------------------------------------

sss0 = np.zeros((3, 62), dtype=complex)
sss5 = np.zeros((3, 62), dtype=complex)
for n1 in range(168):
    for n2 in range(3):
        sss0[n2,:], sss5[n2,:] = seq_gen.sss_gen(n1, n2)
        
end_time = time.time()
#------------------------------------------------------------------

#Output
print("Execution time for pss_gen:", end_time - start_time)   


