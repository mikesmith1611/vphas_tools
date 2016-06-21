from ipyparallel import Client
import os
from astropy.table import Table
import glob
import numpy as np

"""
Calculates shifts for individual pointings
"""


def calibrate(field):
    cmd = 'python /car-data/msmith/tools/calib/calibration.py field={0} apply=True'
    os.system(cmd.format(field))

data = Table.read('bin/vphas_pos_status.fits')

mask = data['finished']
fields = data['Field'][mask]
fields_s = [i.split('_')[1] for i in fields]


def check_done_fields():
    done = glob.glob('tables/*apass-shifts.npy')
    done = [i.split('/')[-1].split('-')[0] for i in done]
    mask = np.in1d(fields, done, invert=True)
    not_done = np.array(fields)[mask]
    woops = []
    missingID = []
    missingSet = []
    for n in not_done:
        d = glob.glob('/car-data/egaps/greimel/VPHAS/merge/*{0}*'.format(n))
        if len(d) >= 2:
            x = [a.endswith('blu.fits') for a in d]
            y = [a.endswith('red.fits') for a in d]
            if not np.array(x).all() and not np.array(y).all():
                woops.append(n)
        else:
            if len(d) == 1:
                if d[0].endswith('red.fits'):
                    missingSet.append('blu')
                elif d[0].endswith('blu.fits'):
                    missingSet.append('red')

            elif len(d) > 1:
                x = [a.endswith('blu.fits') for a in d]
                if np.array(x).all():
                    missingSet.append('red')
                else:
                    missingSet.append('blu')
            else:
                missingSet.append('both')

            missingID.append(n)


def check_done_concats():
    concats = Table.read('bin/vphas_infoblock_p97.fits')
    m = np.in1d(concats['Field'], fields)
    concats = concats[m]
    groupIDs = np.unique(concats['GroupID'])
    done = glob.glob('tables/*apass-shifts.npy')
    done = [i.split('/')[-1].split('-')[0] for i in done]
    badIDs = []
    badFields = []
    goodFields = []
    for ID in groupIDs:
        m1 = concats['GroupID'] == ID
        m2 = np.in1d(concats['Field'][m1], done)
        if len(m2[m2]) != len(m1[m1]):
            badIDs.append(ID)
            badFields.append(list(concats['Field'][m1][-m2]))
            goodFields.append(list(concats['Field'][m1][m2]))
        else:
            print 'huh'
    bad = [badIDs, badFields, goodFields]
    return bad

mycluster = "/home/msmith/.ipython/profile_mpi/security/ipcontroller" \
            "-client.json"

c = Client(mycluster)
v = c[:]

print 'Importing modules...'

with v.sync_imports():
    import os

print "Calibrating..."
sim1 = v.map(calibrate, fields_s, block=True)
