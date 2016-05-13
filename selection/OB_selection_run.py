from ipyparallel import Client
import os
from astropy.table import Table
import glob
import numpy as np


def calibrate(field):
    os.system('python /car-data/msmith/tools/selection/OB_selection_run.py ' +
              field)

data = Table.read('bin/vphas_pos_status.fits')

mask = (data['GAL_LONG'] < 300.) & (data['GAL_LONG'] > 200.) & data['finished']
fields = data['Field'][mask]
fields = [i.split('_')[1] for i in fields]


def check_done():
    done = glob.glob('tables/*apass-shifts.npy')
    done = [i.split('/')[-1].split('-')[0].split('_')[-1] for i in done]
    mask = np.in1d(fields, done, invert=True)
    return np.array(fields)[mask]

mycluster = "/home/msmith/.ipython/profile_mpi/security/ipcontroller" \
            "-calib-client.json"

c = Client(mycluster)
v = c[:]

print 'Importing modules...'

with v.sync_imports():
    import os

print "Calibrating..."
sim1 = v.map(calibrate, fields, block=True)
