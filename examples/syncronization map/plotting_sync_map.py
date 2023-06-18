import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 600
import matplotlib.pylab as pylab

Tt = 100e-9
dtt = 1e-11
# 10 ns at a small time steo

time = np.linspace(0, Tt, int(Tt/dtt))

params = {'legend.fontsize': 'xx-large',
          'figure.figsize': (8.25, 6),
         'axes.labelsize': 'xx-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)

# Will later average from 70ns to 100ns, so need to find where the index of 70ns is in time step
where70ns = np.where(time > 70e-9)[0][0]


# Creating our natural frequency range
natfreqs1 = np.linspace(96e9, 104e9, 300)
natfreqs2 = np.linspace(96e9, 104e9, 300)

# Initializing the frequency difference matrix
freqdiff = np.empty(shape=(len(natfreqs1), len(natfreqs2)))

    # Load freqdiff file
freqdiff = np.genfromtxt('matrix_output.txt')


# Create a threshold matrix at 0.0012 GHz
sync_true = (freqdiff < 0.0012)


fig, ax = plt.subplots(1,1)

c = ax.pcolormesh(natfreqs1 / 1e9, natfreqs2 / 1e9, sync_true, cmap ='Reds_r', vmin=0, vmax=1)
ax.set_xlabel('I1 natural frequency (GHz)')
ax.set_ylabel('I2 natural frequency (GHz)')
ax.text(100, 99.9, "Sync.", color='k', fontsize=18, horizontalalignment='center')
#ax.text(97, 103, "Desync. region", color='w', fontsize=18, horizontalalignment='center')
ax.set_aspect(0.8)

#fig.savefig('synchronization_map_99p5x100p5GHz_range_0p97GHz_coupling')

fig, ax = plt.subplots(1,1)

print(freqdiff)

c = ax.pcolormesh(natfreqs1 / 1e9, natfreqs2 / 1e9, freqdiff, cmap ='Reds', vmin = 0, vmax = np.amax(freqdiff))
ax.set_xlabel('I1 natural frequency (GHz)')
ax.set_ylabel('I2 natural frequency (GHz)')
plt.tight_layout()
cbar = plt.colorbar(c, pad = 0.04)
cbar.set_label('Frequency difference (GHz)', rotation=270, labelpad=+27)
ax.set_aspect(0.8)

#fig.savefig('frequency_map_99p5x100p5GHz_range_0p97GHz_coupling')
