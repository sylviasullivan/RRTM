import numpy as np

base_value = 0.001
factor     = 1.2
lenn       = 32

new_value    = np.zeros((lenn,))
new_value[0] = base_value*factor
for i in np.arange(1,lenn):
    new_value[i] = new_value[i-1]*factor

np.savetxt('IWP_iterator.txt',new_value,fmt='%1.6f')
