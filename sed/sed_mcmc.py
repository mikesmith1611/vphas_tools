import matplotlib
matplotlib.use('Agg')
from astropy.table import Table, join
from astropy.io import fits
import numpy as np
import emcee
import time
from scipy import interpolate
import os
import warnings
import matplotlib.pyplot as plt
import triangle
import vphas

warnings.filterwarnings('ignore')

root = vphas.root


def init(f):

    global field
    field = f
    global star_data
    global star_file
    global path
    global path2
    global Alambda_loc
    global model_loc
    global filterloc
    global filterwavelength

    d = os.getcwd()
    if d.startswith('/local/'):

        star_file = 'outdata/selection/tables/OB_selection_{0}_2mass.fits'
        path = 'outdata/sed/tables/{0}/'
        path2 = 'outdata/sed/plots/{0}/'

        Alambda_loc = root + 'vphas/reddening/Alambda/'
        model_loc = root + 'vphas/reddening/tracks/padova/'
        filterloc = root + 'vphas/reddening/filters/'

    else:
        star_file = ('/car-data/msmith/tools/selection/tables/' +
                     'OB_selection_{0}_2mass.fits')
        Alambda_loc = '../models/reddening/'
        model_loc = '../models/padova/'
        filterloc = '../models/filters/'

        path = 'tables/{0}/'
        path2 = 'plots/{0}/'

    star_data = Table.read(star_file.format(field))
    path = path.format(field)
    path2 = path2.format(field)

    filterwavelength = []
    for i in mod_mag_keys:
        wavelength_lambda, transmissionrate_lambda = np.loadtxt(filterloc + i +
                                                                'filter.txt',
                                                                dtype=float,
                                                                usecols=(0, 1),
                                                                unpack=True)
        filterwavelength.append(np.average(wavelength_lambda,
                                           weights=transmissionrate_lambda,
                                           returned=False))

    global coefs
    coefs = poly_fit(order)

    global extinction_funcs

    extinction_funcs = import_extinction()

sys_u = 0.04
sys_g = 0.03
sys_r = 0.03
sys_i = 0.03
sys_J = 0.03
sys_H = 0.02
sys_Ks = 0.02
sys_errors = [sys_u, sys_g, sys_r, sys_i, sys_J, sys_H, sys_Ks]

ST = {'O2.5': 43651, 'O3.0': 43651, 'O3.5': 42857, 'O4.0': 42857,
      'O4.5': 40862, 'O5.0': 40862, 'O5.5': 39865,
      'O6.0': 38867, 'O6.5': 37870, 'O7.0': 36872, 'O7.5': 35874,
      'O8.0': 34877, 'O8.5': 33879, 'O9.0': 32882, 'O9.5': 31884,
      'O10.0': 30000, 'B1': 25400, 'B2': 22000, 'B3': 18700}

Teff = np.arange(15000, 50200, 200)
RV = np.arange(2.1, 5.41, 0.01)
A0 = np.arange(0, 15.01, 0.01)
dm = np.arange(0, 20, 0.1)

vphas_sys = np.arange(-0.1, 0.1, 0.01)
logTeff = np.log10(Teff)

params = [logTeff, A0, RV, dm]

order = 2.

dat_mag_keys = ['u', 'g', 'r2', 'i', 'Jmag', 'Hmag', 'Kmag']
mod_mag_keys = ['u', 'g', 'r', 'i', 'J', 'H', 'Ks']
fit_mag_keys = ['u', 'g', 'r', 'i', 'J', 'H', 'Ks']
# fit_mag_keys = ['u', 'g', 'r', 'i', 'j_sofi', 'h_sofi', 'k_sofi']
ext_mag_keys = ['Au', 'Ag', 'Ar', 'Ai', 'AJ', 'AH', 'AKs']
err_keys = ['err_u', 'err_g', 'err_r2', 'err_i', 'e_Jmag', 'e_Hmag', 'e_Kmag']
sys_errors = dict(zip(err_keys, sys_errors))
param_keys = ['logteff', 'a0', 'rv', 'dm']
param_labels = [r'$log(Teff) [K]$', r'$A_0$', r'$R_V$', r'$\mu$']
prior_params = dict(zip(param_keys, params))
# coefs = model_poly(mod_mag_keys)
nir_mod = 'Dual'

burn = 200
ndim, nwalkers, nruns = 4, 100, 1000


class sed_model:

    def __init__(self, parameters, mag_keys):

        self.params = parameters
        self.mag_keys = mag_keys
        self.mags = self.model_mags()

    def model_mags(self):

        model_init = teff_func(self.params['logteff'])

        Ax = extinction_sed(self.params['rv'], self.params['a0'])

        model = model_init + Ax + self.params['dm']

        return model


def poly_fit(order):
    filename = (model_loc + 'vst_2mass_sofi_ms.fits')
    data = Table.read(filename)
    coefs = []
    for i, mag in enumerate(fit_mag_keys):
        coefs.append(np.polyfit(data['logTe'], data[mag], order))
    return coefs


def teff_func(logteff):
    sed = [np.polyval(coefs[i], logteff) for i in range(len(mod_mag_keys))]

    return np.array(sed)


class sed_data:

    def __init__(self, data, ID, mag_keys, error_keys):
        mask = (data['ID'] == ID)
        self.data = data[mask]
        self.mag_keys = mag_keys
        self.error_keys = error_keys
        self.mags = self.get_mags()
        self.random_errors = self.get_random_errors()
        self.var = self.get_var()
        self.icov = self.get_icov()
        self.ID = ID

    def get_mags(self):
        mags = np.array([self.data[key][0] for key in self.mag_keys])

        return mags

    def get_random_errors(self):
        errors = np.array([self.data[key][0] for key in self.error_keys])

        return errors

    def get_var(self):
        var_sys = np.array([sys_errors[key] for key in self.error_keys])
        var = self.random_errors ** 2 + np.array(var_sys) ** 2

        return var

    def get_icov(self):
        # cov = np.diag(self.var)

        # icov = np.linalg.inv(cov)

        return 't'  # return icov


class mcmc_results:

    def __init__(self, sampler, param_keys):
        self.param_keys = param_keys
        self.sampler = sampler
        self.flatchain = self.get_flatchains()

    def get_flatchains(self):
        results = [self.sampler.flatchain[:, i]
                   for i in range(len(param_keys))]
        return dict(zip(self.param_keys, results))

    def chain(self, walker):
        chains = [self.sampler.chain[walker][:, i]
                  for i in range(len(param_keys))]
        return dict(zip(self.param_keys, chains))


def uniform_lnprior(x, prior_params):
    for i in x:
        if not (prior_params[i].min() <= x[i] <= prior_params[i].max()):
            return -np.inf
    return 0.0


def distance_lnprior(x, prior_params):
    for i in x:
        if not (prior_params[i].min() <= x[i] <= prior_params[i].max()):
            return -np.inf
    return (10 ** (1 + x['dm'] / 5)) ** 3


def lnprob(x, star_mags, var, mod_mag_keys, prior_params):
    """
    x = [Teff,A0,RV,d]
    star = array of star magnitude values in *mags*
    icov = inverse covariance matrix
    mags = list of magnitudes ['u','g'..etc]

    """
    x = dict(zip(param_keys, x))
    prior = uniform_lnprior(x, prior_params)

    if prior != -np.inf:
        model = sed_model(x, mod_mag_keys).mags
        diff = model - star_mags

        return -0.5 * (np.sum((diff) ** 2 / var)) + prior

    return -np.inf

"""

def A_func():
    R_range = np.arange(2.1, 5.1, 0.01)
    axs = []
    for x in ext_mag_keys:
        ax = []
        for R in R_range:
            A_array = np.copy(Table.read(Alambda_loc
                                         + 'old2/Ax_' + str(R) + '.fits', 1))
            ax.append(A_array[x])
        axs.append(np.array(ax))
    a0s = A_array['A0']

    funcs = []
    for mag in range(len(axs)):
        funcs.append(interpolate.RectBivariateSpline(R_range, a0s, axs[mag]))
    funcs = dict(zip(ext_mag_keys, funcs))

    return funcs


def R_dict():
    R_range = np.arange(2.1, 5.1, 0.01)
    A_array = []
    for R in R_range:
        A_array.append(np.copy(Table.read(Alambda_loc
                                          + 'old2/Ax_' + str(R) + '.fits',
                                          1)))

    Rv = ["%.2f" % x for x in R_range]
    R_dict = dict(zip(Rv, A_array))

    return R_dict


def reddening(rv):
    A = R_dict["%.2f" % rv]

    return A
"""


def import_extinction():

    funcs = []
    for i, j in enumerate(ext_mag_keys):

        hdu = fits.open(Alambda_loc + j + '.fits')[0]

        RV_range = np.arange(hdu.header['RVMIN'],
                             hdu.header['RVMAX'] + 0.01,
                             hdu.header['RVSTEP'])

        A0_range = np.arange(hdu.header['A0MIN'],
                             hdu.header['A0MAX'] + 0.01,
                             hdu.header['A0STEP'])

        func = interpolate.RectBivariateSpline(RV_range, A0_range, hdu.data)

        funcs.append(func)

    funcs = dict(zip(ext_mag_keys, funcs))

    return funcs


def extinction_sed(rv, a0):

    return np.array([float(extinction_funcs[j](rv, a0)) for j in ext_mag_keys])


def mcmcrun(ID, fixed_logteff=True):
    if fixed_logteff:
        mask = star_data['ID'] == ID
        fw_t = star_data['fw_Teff'][mask]
        fw_teff_84 = star_data['fw_Teff_84'][mask]
        fw_teff_16 = star_data['fw_Teff_16'][mask]
        tl_teff = star_data['spec_Teff'][mask]
        tl_teff_84 = tl_teff + star_data['spec_Teff_errup'][mask]
        tl_teff_16 = tl_teff - star_data['spec_Teff_errlo'][mask]

        if fw_t == fw_t:
            Teff = np.linspace(fw_teff_16, fw_teff_84, 20)
            logTeff = np.log10(Teff)
        else:
            Teff = np.linspace(tl_teff_16, tl_teff_84, 20)
            logTeff = np.log10(Teff)
        params = [logTeff, A0, RV, dm]
        prior_params = dict(zip(param_keys, params))

    star = sed_data(star_data, ID, dat_mag_keys, err_keys)

    star_mags = star.mags

    star_var = star.var

    p0 = [np.array([np.random.choice(prior_params[key])
                    for key in param_keys])
          for i in range(nwalkers)]

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob,
                                    args=[star_mags, star_var,
                                          mod_mag_keys, prior_params])

    pos, prob, state = sampler.run_mcmc(p0, burn)

    sampler.reset()

    start = time.time()

    sampler.run_mcmc(pos, nruns)

    end = time.time()

    print end - start

    flatchains = mcmc_results(sampler, param_keys).flatchain

    outdata = Table(flatchains.values(), names=flatchains.keys())

    filename = (path + str(nwalkers) + '_w_' + str(nruns) + '_r_' +
                str(star.ID).strip() + '.fits')

    Table.write(outdata, filename, overwrite=True)

    chains = mcmc_results(sampler, param_keys).chain(nwalkers / 2.)

    outdata = Table(chains.values(), names=chains.keys())

    filename = (path + 'walk_' + str(nwalkers) + '_w_' + str(nruns) + '_r_' +
                str(star.ID).strip() + '.fits')

    Table.write(outdata, filename, overwrite=True)


def mcmctriangle(ID):
    filename = (path + str(nwalkers) + '_w_' + str(nruns) +
                '_r_' + str(ID).strip() + '.fits')

    data = Table.read(filename)

    data_t = np.array([data[key]
                       for key in param_keys]).transpose()

    truths = [np.median(data[key]) for key in param_keys]

    triangle.corner(data_t,
                    labels=param_labels,
                    quantiles=[0.16, 0.5, 0.84], truths=truths)

    filename = (path2 + 'mcmc_' + str(nwalkers) + '_w_' + str(nruns) + '_r_' +
                str(ID).strip() + '.pdf')

    plt.savefig(filename, bbox_inches='tight')
    plt.close()


def bestfits(IDS):
    means = []
    medians = []
    modes = []
    chis = []
    q_84s = []
    q_16s = []

    for i in range(len(IDS)):

        ID = IDS[i]

        filename = (path + str(nwalkers) + '_w_' + str(nruns) +
                    '_r_' + str(ID).strip() + '.fits')

        data = Table.read(filename)

        mean = [np.mean(data[key]) for key in param_keys]
        means.append(mean)

        median = [np.median(data[key]) for key in param_keys]
        medians.append(median)

        def binned_mode(key, bins):
            counts, edges = np.histogram(data[key], bins)

            centers = (edges[:-1] + edges[1:]) / 2.

            mode = centers[counts.argmax()]

            return mode

        mode = [binned_mode(key, 100) for key in param_keys]
        # mode = [np.nan for key in param_keys]
        modes.append(mode)

        q_84 = [np.percentile(data[key], 84) for key in param_keys]
        q_84s.append(q_84)

        q_16 = [np.percentile(data[key], 16) for key in param_keys]
        q_16s.append(q_16)
        averages = [mean, median, mode]

        chis_idv = []

        for j in range(len(averages)):
            star = sed_data(star_data, ID, dat_mag_keys, err_keys)

            star_mags = star.mags

            star_var = star.var

            params = dict(zip(param_keys, averages[j]))

            model = sed_model(params, mod_mag_keys)

            model_mags = model.mags

            chi = np.sum((star_mags - model_mags) ** 2 / star_var)

            chis_idv.append(chi)

        chis.append(chis_idv)

    q_84s = np.array(q_84s)
    q_16s = np.array(q_16s)
    means = np.array(means)
    medians = np.array(medians)
    modes = np.array(modes)
    chis = np.array(chis)

    outdata = [IDS, means[:, 0], means[:, 1], means[:, 2], means[:, 3],
               medians[:, 0], medians[:, 1], medians[:, 2], medians[:, 3],
               modes[:, 0], modes[:, 1], modes[:, 2], modes[:, 3],
               chis[:, 0], chis[:, 1], chis[:, 2], q_16s[:, 0], q_16s[:, 1],
               q_16s[:, 2], q_16s[:, 3], q_84s[:, 0], q_84s[:, 1],
               q_84s[:, 2], q_84s[:, 3]]

    colnames = ['ID', 'logteff_mean', 'a0_mean', 'rv_mean', 'dm_mean',
                'logteff_median', 'a0_median', 'rv_median', 'dm_median',
                'logteff_mode', 'a0_mode', 'rv_mode', 'dm_mode',
                'chi_mean', 'chi_median', 'chi_mode', 'logteff_q_16',
                'a0_q_16', 'rv_q_16', 'dm_q_16', 'logteff_q_84',
                'a0_q_84', 'rv_q_84', 'dm_q_84']

    t1 = Table(outdata, names=colnames)

    table = join(star_data, t1, keys=['ID'])

    filename = path + field + '_averages.fits'

    Table.write(table, filename, overwrite=True)


if __name__ == "__main__":

    from IPython.parallel import Client
    import sys
    import sed_mcmc
    import os
    os.chdir('/car-data/msmith/tools/sed/')
    field = sys.argv[1]
    sed_mcmc.init(field)

    if os.uname()[1] == 'uhppc39.herts.ac.uk':
        c = Client()
        v = c[:]
        print 'Importing Module...'
        v.execute('import sed_mcmc')
        print "Initialising module..."
        v.execute('sed_mcmc.init(' + field + ')')

    else:

        mycluster = "/home/msmith/.ipython/profile_mpi/security/ipcontroller" \
                    "-sedmcmc_" + field + "-client.json"

        c = Client(mycluster)
        v = c[:]
        print 'Importing Module...'
        v.execute('import sed_mcmc')
        print "Initialising module..."
        v.execute('sed_mcmc.init(' + field + ')')

    os.system('mkdir ' + sed_mcmc.path)
    os.system('mkdir ' + sed_mcmc.path2)

    t1 = time.time()

    IDS = list(sed_mcmc.star_data['ID'])

    sim1 = v.map(sed_mcmc.init, [field] * len(v), block=True)
    print "Running MCMC simulation..."

    sim2 = v.map(sed_mcmc.mcmcrun, IDS, block=True)

    print 'Plotting...'

    sim3 = v.map(sed_mcmc.mcmctriangle, IDS, block=True)

    print 'Making Table...'

    sed_mcmc.bestfits(IDS)

    print 'Finished'

    print time.time() - t1
