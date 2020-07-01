import numpy as np

#f = open("tropical_profile_ellingson_250m_formatted.txt", "r")
basedir = '/home/sylvia/Documents/rrtm/'
A = np.genfromtxt(basedir + "output/tropical_profile_ellingson_250m_formatted.txt",dtype=np.float64)
print(A[70,0])
print(A.shape)

Aflip = np.flip(A,axis=0)
print(Aflip[10,0])
print(Aflip.shape)

np.savetxt(basedir + "output/tropical_profile_ellingson_250m_formatted_top2bottom.txt",Aflip)
