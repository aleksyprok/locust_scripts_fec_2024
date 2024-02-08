"""
This module will read in a profiles.CDF and jetto.note file from a JETTO run
and produce the files LOCUST needs to run a simulation.

This includes:
    - a file containing the temperature profile called profile_Ti.dat
      and profile_Te.dat
    - a file containing the density profile called profile_ne.dat
    - a file containing the ion fractions and the alpha power called
      pdep_fi.dat
"""

import sys
import time
import netCDF4  # http://code.google.com/p/netcdf4-python/
import numpy as np
import matplotlib.pyplot as plt

start_time = time.time()

SPR_STRING = 'SPR-045-14'
# MASTER_DIR = (
#     '/home/qt4627/LOCUST/locust_runs/' + SPR_STRING +
#     '/STEP_input_files/' + SPR_STRING
# )
MASTER_DIR = (
    '/home/qt4627/LOCUST/locust_runs/high_vs_low_pedestal/'
    'STEP_input_files/high_vs_low_pedestal'
)
PROFILES_DIR = MASTER_DIR + '/Profiles_' + SPR_STRING

profile_cdf = Dataset(MASTER_DIR + '/Profiles_' + SPR_STRING +
                      '/profiles_' + SPR_STRING + '.CDF',
                      'r')
psin = profile_cdf.variables['XPSI'][-1, :]
ti = profile_cdf.variables['TI'][-1, :]
te = profile_cdf.variables['TE'][-1, :]
ne = profile_cdf.variables['NE'][-1, :]

TI_FILE = PROFILES_DIR + '/profile_Ti.dat'
np.savetxt(TI_FILE, np.transpose([psin, ti]), fmt='%13.5e',
           header=str(len(psin)), comments='')
TE_FILE = PROFILES_DIR + '/profile_Te.dat'
np.savetxt(TE_FILE, np.transpose([psin, te]), fmt='%13.5e',
           header=str(len(psin)), comments='')
NE_FILE = PROFILES_DIR + '/profile_ne.dat'
np.savetxt(NE_FILE, np.transpose([psin, ne]), fmt='%13.5e',
           header=str(len(psin)), comments='')

# Read in the jetto.note file and find the line containing the fusion power
with open(PROFILES_DIR + '/jetto.note', 'r', encoding='utf-8') as f:
    for line in f:
        if 'Fusion Power' in line:
            fusion_power = line.split('=')[-1]
            fusion_power = float(fusion_power[:-3])
            break
print(f'Fusion power: {fusion_power:.2e} GW')
alpha_power = fusion_power * 1e9 * 3.5 / 17.6

NUM_OF_IMPURITIES = 0
nim = []
while True:
    try:
        nim.append(profile_cdf.variables['NIM' +
                                         str(NUM_OF_IMPURITIES + 1)][-1, :])
        NUM_OF_IMPURITIES += 1
    except KeyError:
        break
# Determine impurity atomic numbers and names
impurity_names = []
zia = []
for i in range(NUM_OF_IMPURITIES):
    ziai = profile_cdf.variables[f'ZIA{i+1:d}'][-1, 0]
    ziai = float(ziai.data)
    zia.append(ziai)
    if int(np.round(ziai)) == 2:
        impurity_names.append('AHe4')
    elif int(np.round(ziai)) == 18:
        impurity_names.append('AAr')
    elif 50 <= ziai <= 54:
        impurity_names.append('AXe')
    else:
        print('Error: Unknown impurity')
        print(f'Atomic number: {zia[i]}')
        sys.exit()

print(f'Number of impurities: {NUM_OF_IMPURITIES}')
print(f'Impurity names: {impurity_names}')
print(f'Impurity atomic numbers: {zia}')

# Calculate the ion fractions
ne = profile_cdf.variables['NE'][-1, :]
nid = profile_cdf.variables['NID'][-1, :]
nit = profile_cdf.variables['NIT'][-1, :]
nim = []
for i in range(1, NUM_OF_IMPURITIES + 1):
    nim.append(profile_cdf.variables['NIM' + str(i)][-1, :])

fid = nid[0] / ne[0]
fit = nit[0] / ne[0]
fim = []
for nimi in nim:
    fim.append(nimi[0] / ne[0])

print(f'fid: {fid:.2e}')
print(f'fit: {fit:.2e}')
for i, fimi in enumerate(fim):
    print(f'fim{i + 1}: {fimi:.2e}')

# Check charge neutrality
charge_neutrality = fid + fit + sum(fimi * ziai for fimi, ziai in zip(fim, zia))
print(f'Charge neutrality: {charge_neutrality:.2e}')

# Write the ion fractions, ion names and alpha power to a file
# in the following format:
# number_of_impurities + 2
# alpha_power
# zia fi ion_name
# zia fi ion_name
# ...
PDEP_FI_FILE = PROFILES_DIR + '/pdep_fi.dat'
with open(PDEP_FI_FILE, 'w', encoding='utf-8') as f:
    f.write(f'{NUM_OF_IMPURITIES + 2:d}\n')
    f.write(f'{alpha_power:.5e}\n')
    f.write(f'{1:11.8f} {fid:11.8f} AD\n')
    f.write(f'{1:11.8f} {fit:11.8f} AT\n')
    for i, (ziai, fimi) in enumerate(zip(zia, fim)):
        f.write(f'{ziai:11.8f} {fimi:11.8f} {impurity_names[i]}\n')

end_time = time.time()
print(f"Time taken: {end_time - start_time:.2e} seconds")

plt.show()
