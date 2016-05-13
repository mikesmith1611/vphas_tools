import os
import glob
import subprocess
from astropy.table import Table

data = Table.read('bin/carina_fields.fits')
fields = data['Field']


def check_completion():
    files = glob.glob('/home/msmith/sedmcmc_????.o*')

    todo = []
    for f in files:

        line = subprocess.check_output(['tail', '-1', f])

        if line == 'Terminated\n':

            field = f.split('.')[0].split('_')[1]
            todo.append(field)

    return todo


def submit_job(fields):

    for f in fields:

        print f

        os.system('qsub /car-data/msmith/carina/outdata/selection/tables/' +
                  f + '_mcmc.pbs')

        old_logs = glob.glob('/home/msmith/sedmcmc_' + f + '.*')

        for log in old_logs:

            os.system('rm ' + log)


if __name__ == "__main__":

    fields = [f.split('_')[1] for f in fields]
    submit_job(fields)
