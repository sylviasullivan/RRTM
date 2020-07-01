import matplotlib.pyplot as plt
import numpy as np

# Pressure at half levels
basedir = '/home/sylvia/Documents/rrtm/'
ellingson = np.genfromtxt(basedir + 'output/tropical_profile_ellingson_250m_formatted_top2bottom.txt')
pp_fl = ellingson[:,1]
pp_hl = np.zeros((82,))
for i in np.arange(1,81):
    pp_hl[i] = (pp_fl[i-1] + pp_fl[i])/2.0
pp_hl[0] = pp_fl[0] - (pp_hl[1] - pp_fl[0])
pp_hl[81] = pp_fl[80] + (pp_fl[80] - pp_hl[80])

# Heat capacity [J kg-1 K-1]
cp = 1.08*10**(3)
# Gravity [m s-2]
g = 9.8
# Time of day at which clr-sky profile was calculated
tt = np.asarray(['01','02','03','04','05','06','07','08','09',
       '10','11','12','13','14','15','16','17','18','19','20',
       '21','22','23'])

#farbe = ['red','orange','gold','green','blue','purple']
#lbl = ['261 K','262 K','264 K','266 K','267 K','269 K']

fs = 13
fig = plt.figure()
for c,tk in enumerate(tt):
    lw = np.genfromtxt(basedir + 'lwflxatm-test_diurnal' + tk + '.txt')
    # First column is all-sky. Second column is clear-sky. 
    swclr = lw[1]

    # Calculate the longwave heating rate. 
    # Factor of 3600 converts from K s-1 to K day-1.
    #H = g/cp*np.gradient(lwcld,pp_hl*100.)*3600
    H = g/cp*np.gradient(swclr,pp_hl*100.)*3600
    plt.plot(H,pp_hl,linewidth=1.25) #color=farbe[c],linewidth=1.25,label=lbl[c])

    plt.ylabel('Pressure [hPa]',fontsize=fs)
    plt.xlabel(r'Clr-sky LW heating rate [K day$^{-1}$]',fontsize=fs)

plt.legend()
plt.gca().invert_yaxis()
#fig.savefig('lwheating_clrsky.pdf',bbox_inches='tight')
plt.show()
