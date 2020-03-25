import matplotlib.pyplot as plt
import numpy as np
import sys

namefile = sys.argv[1]
print("file name: ", namefile)
loaded = np.loadtxt(namefile, skiprows = 1, usecols = range(2,5) )

measured = loaded[:,0]
modelled = loaded[:,1]
residual = loaded[:,2]
numparam = residual.shape[0]
onesvect = np.linspace(1,numparam,numparam)
maxD = max(max(measured), max(modelled))
minD = min(min(measured), min(modelled))

font = {'family':'serif','style':'normal','weight':'normal', 'size':12}
plt.rc('font', **font)


plt.figure(1, figsize = (12,8))
plt.subplot(121)
plt.plot(measured[:306],modelled[:306],'og')
plt.xlabel('LabData [ohm]', fontsize = 14)
plt.ylabel('SynData [ohm]', fontsize = 14)
plt.xlim(minD-2, maxD+2)
plt.ylim(minD-2, maxD+2)
ax= plt.gca()
ax.set_aspect('equal', adjustable='box')
ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")

plt.subplot(122)
plt.plot(onesvect,measured,'or', label = 'Lab [ohm]')
plt.plot(onesvect,modelled,'ob', label = 'Syn [ohm]')
plt.xlabel('Data number', fontsize = 14)
plt.ylabel('Resistance [ohm]', fontsize = 14)
plt.legend()
plt.savefig('CheckPEST.png')
plt.show()
