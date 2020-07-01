import numpy as np
import matplotlib.pyplot as plt
import sys
exec(open("qvsat.py").read())

ellingson = np.genfromtxt('tropical_profile_ellingson_250m_formatted_top2bottom.txt')
z = ellingson[:,0]
T = ellingson[:,2]
P = ellingson[:,1]
qv = ellingson[:,5]
RH = qv/qvsat(T)

# Calculate lapse rates from above.
dTdz = np.gradient(T,z)
# Define the tropopause by the cold-point temperature.
tropopause = np.argmin(T)
MLR = -1.*np.nanmean(dTdz[tropopause:])
# Calculate the mean RH in the troposphere.
RHmean = np.nanmean(RH[tropopause:])

fig,ax = plt.subplots(nrows=1,ncols=4,figsize=[12,4])
fs = 13


ax[0].plot(T,z,linewidth=1.25)
ax[0].text(0.05,0.8,'(a)',fontsize=fs+2,weight='bold',transform=ax[0].transAxes)
ax[0].text(0.05,0.35,r'$T_s$ = {} K'.format(T[-1]),transform=ax[0].transAxes)
ax[0].text(0.05,0.25,r'$<\Gamma >$ = {:.2f}'.format(MLR) + r' K km$^{-1}$',transform=ax[0].transAxes)
ax[0].set_ylabel('Altitude [km]',fontsize=fs-1)
ax[0].set_xlabel('Temperature [K]',fontsize=fs-1)
ax[0].set_ylim([0,15])

ax[1].plot(P,z)
ax[1].text(0.05,0.8,'(b)',fontsize=fs+2,weight='bold',transform=ax[1].transAxes)
ax[1].set_xlabel('Pressure [hPa]',fontsize=fs-1)
ax[1].set_ylim([0,15])

ax[2].plot(qv*1000.,z)
ax[2].text(0.05,0.85,'(c)',fontsize=fs+2,weight='bold',transform=ax[2].transAxes)
ax[2].set_xlabel(r'Specific humidity [g kg$^{-1}$]',fontsize=fs-1)
ax[2].set_ylim([0,15])

ax[3].plot(RH,z)
ax[3].text(0.05,0.85,'(d)',fontsize=fs+2,weight='bold',transform=ax[3].transAxes)
ax[3].text(0.05,0.35,r'$<$RH$>$ = '
                     '\n'
                     '{:.3f}'.format(RHmean),transform=ax[3].transAxes)
ax[3].set_xlabel('Relative humidity',fontsize=fs-1)
ax[3].set_ylim([0,15])
ax[3].set_xlim([0,1])

fig.savefig('ellingson_profiles.pdf',bbox_inches='tight')
plt.show()
