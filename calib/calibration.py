import numpy as np
from astropy.table import Table, vstack
import glob
import os
import scipy.optimize as opt
import fnmatch
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.colors as colors
import colormaps as cmaps
from matplotlib import rc, rcParams
import matplotlib.path as pltpath
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
rc('text', usetex=True)
rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{color}"]

# Change all fonts to 'Comptuer Modern'
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})

# warnings.filterwarnings('ignore')

# Location of where to store the merged VPHAS+ data
field_loc = '/car-data/msmith/fields/'
wd = '/car-data/msmith/tools/calib/'
stilts = 'java -jar /home/msmith/bin/topcat/topcat-*.jar -stilts'


def density_plot(data, cx, cy, cxshift, cyshift, ax, bins):

    mx = []
    my = []

    sp = "G0_V"
    for i in cx:
        i = i.split('_')[0]
        if i == 'r2' or i == 'r1':
            mx.append('r')

        elif i == 'g':
            mx.append('g')

        elif i == 'u':
            mx.append('u')

        elif i == 'Ha':
            mx.append('Ha')
            sp = 'A2_V'
        else:
            mx.append(i)
        # if i == 'i':
        #     sp = 'A2_V'
    for i in cy:
        i = i.split('_')[0]
        if i == 'r1' or i == 'r2':
            my.append('r')
            sp = "O6_V"

        elif i == 'g':
            my.append('g')

        elif i == 'u':
            my.append('u')

        elif i == 'Ha':
            my.append('Ha')
            # sp = 'A2_V'
        else:
            my.append(i)

    reddata = Table.read('../models/reddening/Rv_3.1.fits')

    mask = (reddata['spectype'] == sp)
    model_GOV = reddata[mask]

    GOV_gr = model_GOV[mx[0]] - model_GOV[mx[1]]
    GOV_ug = model_GOV[my[0]] - model_GOV[my[1]]

    x = (data[cx[0]] - data[cx[1]]) + cxshift

    y = (data[cy[0]] - data[cy[1]]) + cyshift

    mask = ((x == x) & (y == y))

    x = np.copy(x[mask])
    y = np.copy(y[mask])

    H, xedges, yedges = np.histogram2d(y, x, bins=bins)
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]

    levels = 10 ** np.linspace(0, np.log10(H.max()), 10)

    norm = colors.LogNorm(vmin=levels.min(), vmax=levels.max())
    cf = ax.contourf(H, extent=extent, levels=levels,
                     cmap=cmaps.viridis, norm=norm)

    for a in np.arange(0, 6, 1):
        mask = (reddata['A0'] == a)
        ZAMS = reddata[mask]
        mask = np.argsort(ZAMS[mx[0]] - ZAMS[mx[1]])
        ZAMS = ZAMS[mask]

        ZAMS_gr = ZAMS[mx[0]] - ZAMS[mx[1]]
        ZAMS_ug = ZAMS[my[0]] - ZAMS[my[1]]
        ax.plot(ZAMS_gr, ZAMS_ug, color='k', linewidth=0.5, alpha=0.5)
    ax.plot(GOV_gr, GOV_ug, color='k', linewidth=0.5, linestyle='dashdot',
            alpha=0.8)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", pad=0.05)
    cb = plt.colorbar(cf, cax=cax, orientation="horizontal")

    vals = []

    for i, c in enumerate(levels):
        if i != 0:
            mask1 = H <= levels[i]
            mask2 = H > levels[i - 1]
            mask = mask1 & mask2
        else:
            mask = H <= levels[i]

        vals.append(np.sum(H[mask]))

    labels = []

    for v in vals:
        a, b = "{0:.1e}".format(v).split("e+")
        c = r'$' + a + r'\times 10^{' + b[1] + '}$'
        labels.append(c)

    labels = ["{0}".format(int(x)) for x in levels]
    # cbax = cb.ax.twiny()
    cb.ax.set_xlabel("Density ($N_{Stars}/0.01$ mag$^2$)", fontsize=8,
                     labelpad=-28)

    cb.ax.tick_params(labelsize=4, which="both", pad=0.01,
                      top='off', bottom='off',
                      labeltop=True, labelbottom=False)
    cb.ax.set_xticklabels(labels, fontsize=8)
    #cb.ax.get_xticklabels()[0].set_horizontalalignment('left')
    #cb.ax.get_xticklabels()[-1].set_horizontalalignment('right')
    # cb.ax.set_xticks([], visible=False)
    # cb.ax.set_xticklabels([], visible=False)

    return ax

def getmodelcolours(**kwargs):

    d = {'R_V': 3.1, 'ms': False, 'A0': 0, 'spectype': 'G0_V'}

    for kwarg in kwargs:
        d[kwarg] = kwargs[kwarg]

    reddata = Table.read('../models/reddening/Rv_{0}.fits'.format(d['R_V']))

    if d['ms']:
        mask = (reddata['A0'] == d['A0'])
        ms = reddata[mask]
        return ms
    if not d['ms']:
        mask = (reddata['spectype'] == d['spectype'])
        redvector = reddata[mask]
        return redvector

###############################################################################


def calib_polygon():
    """
    :return:
    """
    G0model = getmodelcolours(spectype=['G0_V'], R_V=3.1)
    mask = (G0model['A0'] < 2.7)
    G0model = G0model[mask]

    ZAMS = getmodelcolours(ms=True)
    mask = np.argsort(ZAMS['g'] - ZAMS['r'])
    ZAMS = ZAMS[mask]
    mask = (ZAMS['g'] - ZAMS['r'] > 0.65) & (ZAMS['g'] - ZAMS['r'] < 1.437)
    ZAMS = ZAMS[mask]

    s1 = zip(G0model['g'] - G0model['r'], G0model['u'] - G0model['g'])
    s2 = zip(ZAMS['g'] - ZAMS['r'], ZAMS['u'] - ZAMS['g'])

    poly = s1 + s2[::-1]

    return pltpath.Path(poly)

bbPath = calib_polygon()

###############################################################################


def crossmatch(in1, in2, out, **kwargs):
    """
    ====================================================================
    crossmatch - joins red and blu VPHAS+ tables
    ====================================================================
    Variables

    in1 - input file location table 1 (string)

    in2 - input file location table 2 (string)

    out - output file location (string)

    **kwargs

    params - cross match radius(arcsec) (default=0.5)

    :param in1:
    :param in2:
    :param out:
    :param kwargs:
    """
    d = {'params': 0.5, 'join': '1or2', 'fixcols': 'dups'}
    for kwarg in kwargs:
        d[kwarg] = kwargs[kwarg]
    os.system('rm ' + out)
    os.system(stilts + ' tmatch2 ' + 'in1=' + str(in1) + ' ' +
              'in2=' + str(in2) +
              ' ' + 'out=' + str(out) +
              ' matcher=sky values1="RA DEC" values2="RA DEC" ' +
              'params=' + str(d['params']) + ' join=' + d['join'] +
              ' fixcols=' + d['fixcols'])

###############################################################################


class vphas_field:
    def __init__(self, field):

        self.field = field
        self.a2v = {'r1': 'rmag', 'r2': 'rmag', 'g': 'gmag', 'i': 'imag'}
        self.lims = {'r1': (13.5, 16.), 'r2': (13.5, 16.),
                     'g': (14, 16.), 'i': (13.5, 15.),
                     'u': (13, 20)}
        self.ab_to_vega = {'r1': -0.136, 'r2': -0.136, 'g': 0.123, 'i': -0.373}
        self.bbPath = bbPath

    def download(self):

        filenames = glob.glob('/car-data/egaps/greimel/VPHAS/merge/*')
        prefix = ''

        x = fnmatch.filter(filenames, '*vphas_' + self.field + '*')

        blu = fnmatch.filter(x, '*blu*')

        red = fnmatch.filter(x, '*red*')

        if len(blu) > 1:

            dates = np.array([int(blu[i].split(self.field)[-1].split('-')[1])
                              for i, j in enumerate(blu)])

            blu_file = prefix + blu[dates.argmax()]
            os.system('rsync -avuhP ' + blu_file + ' ' +
                      field_loc + 'allfields/' +
                      blu_file.split('/')[-1])

        elif len(blu) == 1:

            blu_file = prefix + blu[0]

            os.system('rsync -avuhP ' + blu_file + ' ' +
                      field_loc + 'allfields/' +
                      blu_file.split('/')[-1])

        else:
            print 'file does not exist'

        if len(red) > 1:

            dates = np.array([int(red[i].split(self.field)[-1].split('-')[1])
                              for i, j in enumerate(red)])

            red_file = prefix + red[dates.argmax()]

            os.system('rsync -avuhP ' + red_file + ' ' +
                      field_loc + 'allfields/' +
                      red_file.split('/')[-1])

        elif len(red) == 1:

            red_file = prefix + red[0]
            os.system('rsync -avuhP ' + red_file + ' ' +
                      field_loc + 'allfields/' +
                      red_file.split('/')[-1])
        else:
            print 'file does not exist'

        download(self.field)

        files = glob.glob(field_loc + 'allfields/vphas_' +
                          self.field + '-*.fits')

        print len(files)

        if len(files) == 2:

            blu = fnmatch.filter(files, '*blu*')[0]
            red = fnmatch.filter(files, '*red*')[0]

            blu_data = Table.read(blu)

            rcols = fnmatch.filter(blu_data.keys(), '*r_1*')
            for c in rcols:
                newc = c.replace('r_1', 'r2_1')
                blu_data.rename_column(c, newc)

            rcols = fnmatch.filter(blu_data.keys(), '*r_2*')
            for c in rcols:
                newc = c.replace('r_2', 'r2_2')
                blu_data.rename_column(c, newc)

            blu_data.write(blu, overwrite=True)

            red_data = Table.read(red)

            rcols = fnmatch.filter(red_data.keys(), '*r_1*')

            for c in rcols:
                newc = c.replace('r_1', 'r1_1')
                red_data.rename_column(c, newc)

            rcols = fnmatch.filter(red_data.keys(), '*r_2*')

            for c in rcols:
                newc = c.replace('r_2', 'r1_2')
                red_data.rename_column(c, newc)

            red_data.write(red, overwrite=True)

            file = (field_loc + 'allfieldsxmatched/vphas_' +
                    self.field + '.fits')

            crossmatch(blu, red, file)

            os.system(stilts + ' tpipe in=' + str(file) +
                      ' cmd=' + "'addcol" + ' RA ' +
                      '"NULL_RA_1 ? RA_2:RA_1"' + "'" +
                      ' cmd=' + "'addcol" + ' DEC ' +
                      '"NULL_DEC_1 ? DEC_2:DEC_1"' + "'" +
                      ' out=' + str(file) + ' ofmt=fits')

            for fi in files:
                os.system('rm ' + fi)

        else:
            print 'Missing files!'
            return 'Missing files!'

    def apass_clean(self, data, band, offset):

        lims = self.lims
        a2v = self.a2v

        col = band + '_' + offset
        er = 'err_' + col

        mask1 = (data[col] == data[col]) & (data[a2v[band]] == data[a2v[band]])

        mask2 = ((data[col] > lims[band][0]) & (data[col] < lims[band][1]) &
                 (data[er] < 0.1) & (data[a2v[band]] < 99))

        mask = mask1 & mask2

        return data[mask]

    def plot_apass(self, ax, data, diff, band, shift):

        magstring = band.split('_')[0]

        ax.plot(data[band], diff,
                mfc='#888888', mec='#888888', marker='.',
                linestyle='None')

        ax.hlines(shift, self.lims[magstring][0], self.lims[magstring][1],
                  linestyles='dashed', zorder=4)

        ax.set_xlim(self.lims[magstring][0], self.lims[magstring][1])

        ax.set_title(r'$' + magstring + r'$')

        ax.set_xlabel(r'$' + magstring + r'$')

        ax.set_ylabel(r'$' + magstring +
                      r'(VPHAS+) - ' + magstring + r'(APASS)$')

    def apass_calibrate(self,
                        recalibrate={'r1': False, 'r2': False,
                                     'g': False, 'i': False, 'u': False}):

        if np.any(recalibrate.values()):
            shift_data = Table.read('tables/shifts_all.fits')

        vphas_file = (field_loc + 'allfieldsxmatched/vphas_' +
                      self.field + '.fits')

        outfile = (wd + 'tables/vphas_' +
                   self.field + '-apass.fits')

        if not os.path.exists(outfile):
            crossmatch(vphas_file, '../models/apass/' +
                       'apass-dr9-plane-vphas.fits',
                       outfile, join='1and2')

        data = Table.read(outfile)
        vdata = Table.read(vphas_file)

        shifts = {'r1': [], 'r2': [], 'g': [], 'i': [], 'u': []}
        sigmas = {'r1': [], 'r2': [], 'g': [], 'i': []}
        subplots = {'r1': (1, 0), 'r2': (1, 1), 'g': (2, 0), 'i': (2, 1),
                    'u_g_b': (0, 2), 'u_g': (0, 3), 'r_ha_b': (1, 2),
                    'r_ha': (1, 3), 'r_i_b': (2, 2), 'r_i': (2, 3)}

        for i in ['1', '2']:
            fig = plt.figure(figsize=(16, 12))

            for b in self.a2v.keys():

                d = self.apass_clean(data, b, i)

                delta = d[b + '_' + i] - (d[self.a2v[b]] + self.ab_to_vega[b])

                mask_delta = (np.abs(delta - np.median(delta)) >
                              (2 * np.std(delta)))

                ax = plt.subplot2grid((3, 4), subplots[b])

                if len(mask_delta[mask_delta]) != 0:
                    shift = np.median(delta[-mask_delta])
                    self.plot_apass(ax, d[-mask_delta],
                                    delta[-mask_delta], b + '_' + i, shift)
                else:
                    shift = np.median(delta)
                    self.plot_apass(ax, d,
                                    delta, b + '_' + i, shift)

                if recalibrate[b]:
                    shift = np.median(shift_data[b])
                    lim = ax.get_xlim()
                    ax.hlines(shift, lim[0], lim[1], color='r',
                              linestyles='dashed', zorder=4)

                shifts[b].append(-shift)
                sigmas[b].append(np.std(delta))

            ub = 'u' + '_' + i

            mask1 = vdata['err_' + ub] < 0.1
            mask2 = vdata[ub] > self.lims['u'][0]
            mask3 = vdata[ub] < self.lims['u'][1]

            vdatau = vdata[mask1 & mask2 & mask3]

            g = vdatau['g' + '_' + i] + shifts['g'][int(i) - 1]

            r = vdatau['r2' + '_' + i] + shifts['r2'][int(i) - 1]

            u = vdatau[ub]

            bbPath = self.bbPath

            def density_func(ushift):
                inside = bbPath.contains_points(zip(g - r, u - g + ushift))
                return -len(inside[inside])

            res = opt.minimize_scalar(density_func, tol=0.01)
            shift_best = res.x

            x = g - r

            y = u - g + shift_best

            mask = ((x == x) & (y == y))

            x = np.copy(x[mask])
            y = np.copy(y[mask])
            xlims = [0., 2.0]
            ylims = [2.0, -0.5]
            binsx = np.arange(min(xlims), max(xlims) + 0.01, 0.01)
            binsy = np.arange(min(ylims), max(ylims) + 0.01, 0.01)
            H, yedges, xedges = np.histogram2d(y, x, bins=[binsy, binsx])
            xcent = (xedges[1:] - xedges[:-1]) / 2 + xedges[:-1]
            ycent = (yedges[1:] - yedges[:-1]) / 2 + yedges[:-1]

            xarr, yarr = np.meshgrid(xcent, ycent)

            levels = 10 ** np.linspace(0, np.log10(H.max()), 10)
            mask1 = H >= levels[-2]
            xmax, ymax = [xarr[mask1], yarr[mask1]]
            mask2 = bbPath.contains_points(zip(xmax, ymax))
            n_in = len(mask2[mask2])
            n_out = len(mask2[-mask2])

            if n_out > n_in:
                print 'u band additional shift'
                recalibrate['u'] = True
                model_GOV = getmodelcolours(spectype=['G0_V'], R_V=3.1)
                mask = np.argsort(model_GOV['g'] - model_GOV['r'])
                GOV_gr = model_GOV['g'] - model_GOV['r']
                GOV_ug = model_GOV['u'] - model_GOV['g']
                z = np.polyfit(GOV_gr[mask], GOV_ug[mask], 3)

                def obj_count(unewshift):
                    mask = (ymax + unewshift < np.lib.polyval(z, xmax))
                    return abs(len(mask[mask]) - len(mask[-mask]))

                its = []
                uranges = np.arange(-0.2, 0.21, 0.01)
                for nu in uranges:
                    its.append(obj_count(nu))
                unewshift = uranges[np.argmin(its)]

                shift_best += unewshift

            shifts['u'].append(shift_best)

            for k in ['u_g_b', 'u_g']:

                ax = plt.subplot2grid((3, 4), subplots[k])

                if k.endswith('g'):
                    s1 = shifts['g'][int(i) - 1] - shifts['r2'][int(i) - 1]
                    s2 = shifts['u'][int(i) - 1] - shifts['g'][int(i) - 1]

                else:
                    s1 = 0
                    s2 = 0

                xlims = [0., 2.0]
                ylims = [2.0, -0.5]
                binsx = np.arange(min(xlims), max(xlims) + 0.01, 0.01)
                binsy = np.arange(min(ylims), max(ylims) + 0.01, 0.01)

                ax = density_plot(vdatau, ('g_' + i, 'r2_' + i),
                                  ('u_' + i, 'g_' + i),
                                  s1, s2, ax, [binsy, binsx])

                bb = bbPath.to_polygons()

                p = mpl.patches.Polygon(bb[0], closed=True,
                                        facecolor='None',
                                        edgecolor='k', zorder=5)
                ax.add_artist(p)

                ax.set_ylim(ylims)
                ax.set_xlim(xlims)
                ax.set_xlabel(r'$g-r$', fontsize=16)
                ax.set_ylabel(r'$u-g$', fontsize=16)

            mask1 = vdata['Av_conf_Ha_' + i] > 95.
            mask2 = vdata['err_Ha_' + i] < 0.1

            vdataha = vdata[mask1 & mask2]

            for k in ['r_ha_b', 'r_ha']:
                ax = plt.subplot2grid((3, 4), subplots[k])

                if k.endswith('ha'):
                    s1 = shifts['r1'][int(i) - 1] - shifts['i'][int(i) - 1]
                else:
                    s1 = 0

                xlims = [-0.2, 2.25]
                ylims = [-0.5, 1.2]
                binsx = np.arange(min(xlims), max(xlims) + 0.01, 0.01)
                binsy = np.arange(min(ylims), max(ylims) + 0.01, 0.01)
                ax = density_plot(vdataha,
                                  ['r1_' + i, 'i_' + i],
                                  ['r1_' + i, 'Ha_' + i],
                                  s1, 0, ax, [binsy, binsx])

                ax.set_ylim(ylims)
                ax.set_xlim(xlims)
                ax.set_xlabel(r'$r-i$', fontsize=16)
                ax.set_ylabel(r'$r-ha$', fontsize=16)

            for k in ['r_i_b', 'r_i']:
                ax = plt.subplot2grid((3, 4), subplots[k])

                if k.endswith('i'):
                    s1 = shifts['r1'][int(i) - 1] - shifts['i'][int(i) - 1]
                    s2 = shifts['g'][int(i) - 1] - shifts['r2'][int(i) - 1]
                else:
                    s1 = 0
                    s2 = 0
                ylims = [-0.25, 3.0]
                xlims = [-0.25, 2.0]
                binsx = np.arange(min(xlims), max(xlims) + 0.01, 0.01)
                binsy = np.arange(min(ylims), max(ylims) + 0.01, 0.01)

                ax = density_plot(vdatau,
                                  ['r1_' + i, 'i_' + i],
                                  ['g_' + i, 'r2_' + i],
                                  s1, s2, ax, [binsy, binsx])

                ax.set_ylim(ylims)
                ax.set_xlim(xlims)
                ax.set_xlabel(r'$r-i$', fontsize=16)
                ax.set_ylabel(r'$g-r$', fontsize=16)
            plt.tight_layout()
            txt = []
            cols = ['u', 'g', 'r1', 'r2', 'i']
            for s in cols:
                t = '{0:.2f}'.format(shifts[s][int(i) - 1])
                txt.append(t)

            ax = plt.subplot2grid((3, 4), (0, 0), colspan=2)
            ytable = ax.table(cellText=[txt], colLabels=cols,
                              loc='center', cellLoc='center', fontsize=20)

            for f, s in enumerate(cols):
                if recalibrate[s]:
                    ytable._cells[1, f]._text.set_color('r')
            ytable.scale(1, 4)
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            ax.axis('off')
            plt.subplots_adjust(top=0.92, bottom=0.05)
            ax.text(0.45, 0.8, r'\underline{Shifts}',
                    ha='center', va='center', fontsize=20)
            fig.text(0.25, 0.65, r'\underline{APASS Comparsion}',
                     ha='center', va='center', fontsize=20)
            fig.text(0.65, 0.97, r'\underline{Before Calibration}',
                     ha='center', va='center', fontsize=20)
            fig.text(0.89, 0.97, r'\underline{After Calibration}',
                     ha='center', va='center', fontsize=20)

            plt.savefig(wd + 'plots/vphas_' +
                        self.field + '_' + i + '.png', dpi=100)

            plt.close()

        np.save(outfile[:-5] + '-shifts.npy', np.array([shifts]))
        np.save(outfile[:-5] + '-sigmas.npy', np.array([sigmas]))
        return shifts

    def sort(self):

        sort(field_loc + 'allsortedfields/vphas_' +
             self.field + '.fits',
             field_loc + 'allsortedfields/vphas_' +
             self.field + '.fits')

        print '\n' + 'Sorting narrow band...'

        merge(field_loc + 'allsortedfields/vphas_' +
              self.field + '.fits',
              field_loc + 'allsortedfields/vphas_' +
              self.field + '.fits')

        os.system('mv ' + field_loc + 'allsortedfields/vphas_' +
                  self.field + '.fits ' + wd +
                  'fields/vphas_' +
                  self.field + '.fits')

        cc_plots(self.field)

    def aply_shifts(self):

        data = Table.read(field_loc + 'allfieldsxmatched/vphas_' +
                          self.field + '.fits')

        bands = ['Ha', 'r1', 'r2', 'g', 'i', 'u']

        shifts = Table.read(wd + 'tables/shifts_all.fits')

        for band in bands:

            for o in ['1', '2']:
                mask = (shifts['field'] == self.field + '_' + o)
                shift = shifts[mask][band + '_c'][0]
                data[band + '_' + o] = data[band + '_' + o] + shift

            if band == 'Ha':
                mask = (shifts['field'] == self.field + '_2')
                shift = shifts[mask][band + '_c'][0]
                data[band + '_3'] = data[band + '_3'] + shift

        data.write(field_loc + 'allsortedfields/vphas_' +
                   self.field + '.fits',
                   overwrite=True)

        # os.system('rm ' + field_loc + 'allfieldsxmatched/vphas_' +
        #           self.field + '.fits')

###############################################################################


def calib_polygon2(cx, cy):
    """
    :return:
    """
    G0model = getmodelcolours(spectype=['G0_V'], R_V=3.1)
    mask = (G0model['A0'] < 2.7)
    G0model = G0model[mask]

    ZAMS = getmodelcolours(ms=True)
    mask = np.argsort(ZAMS['g_r'])
    ZAMS = ZAMS[mask]
    mask = (ZAMS[cx] > 0.65) & (ZAMS[cy] < 1.437)
    ZAMS = ZAMS[mask]

    s1 = zip(G0model[cx], G0model[cy])
    s2 = zip(ZAMS[cx], ZAMS[cy])

    poly = s1 + s2[::-1]

    return pltpath.Path(poly)


def cc_plots(field):

    datain = Table.read('fields/vphas_{0}.fits'.format(field))
    shiftsa = Table.read('tables/shifts_all.fits')
    cx = [('g', 'r2'), ('r1', 'i'), ('r1', 'i')]
    cy = [('u', 'g'), ('r1', 'Ha'), ('g', 'r2')]
    xlims = [(0., 2.0), (-0.25, 2.25), (-0.25, 2.0)]
    ylims = [(2.0, -0.5), (-0.5, 1.2), (-0.25, 3.0)]
    splots = [(1, 2), (3, 4), (5, 6)]

    for j in ['1', '2']:
        mask = shiftsa['field'] == field + '_' + j
        shifts = shiftsa[mask]
        mask2 = (datain['Av_conf_Ha_' + j] > 95.) & (datain['err_Ha_' + j] < 0.1)
        mask1 = (datain['u_' + j] > 13) & (datain['u_' + j] < 20) & (datain['err_u_' + j] < 0.1)
        mask3 = (datain['g_' + j] > 13) & (datain['g_' + j] < 20) & (datain['err_g_' + j] < 0.1)
        masks = [mask1, mask2, mask3]
        plt.figure(figsize=(8, 12))
        for i in range(len(cx)):
            data = datain[masks[i]]
            cx2 = (cx[i][0] + '_' + j, cx[i][1] + '_' + j)
            cy2 = (cy[i][0] + '_' + j, cy[i][1] + '_' + j)
            ax1 = plt.subplot(320 + splots[i][0])
            binsx = np.arange(min(xlims[i]), max(xlims[i]) + 0.01, 0.01)
            binsy = np.arange(min(ylims[i]), max(ylims[i]) + 0.01, 0.01)
            ax1 = density_plot(data, cx2, cy2,
                               -float(shifts[cx[i][0] + '_c'] -
                                      shifts[cx[i][1] + '_c']),
                               -float(shifts[cy[i][0] + '_c'] -
                                      shifts[cy[i][1] + '_c']),
                               ax1, [binsy, binsx])
            ax1.set_xlabel('{0} - {1}'.format(cx[i][0], cx[i][1]))
            ax1.set_ylabel('{0} - {1}'.format(cy[i][0], cy[i][1]))
            ax1.set_xlim(xlims[i])
            ax1.set_ylim(ylims[i])
            ax2 = plt.subplot(320 + splots[i][1])
            ax2 = density_plot(data, cx2, cy2,
                               0,
                               0,
                               ax2, [binsy, binsx])
            ax2.set_xlabel('{0} - {1}'.format(cx[i][0], cx[i][1]))
            ax2.set_ylabel('{0} - {1}'.format(cy[i][0], cy[i][1]))
            ax2.set_xlim(xlims[i])
            ax2.set_ylim(ylims[i])
            if i == 0:
                x = bbPath.to_polygons()

                p = mpl.patches.Polygon(x[0], closed=True,
                                        facecolor='None',
                                        edgecolor='k', zorder=5)
                ax1.add_artist(p)
                ax2.add_artist(p)

        plt.tight_layout()
        plt.savefig('plots/vphas_{0}_{1}_concat_calib.png'.format(field, j),
                    bbox_inches='tight')


def make_before_after_calib():

    fields = Table.read('bin/carina_fields.fits')['Field']
    shifts = Table.read('outdata/calib/tables/shifts.fits')
    tables_b = []
    tables_a = []
    for f in fields:
        f = f.split('_')[1]
        print f
        data_a = Table.read('outdata/calib/fields/vphas_{0}.fits'.format(f))

        mask = ((data_a['g'] > 13) & (data_a['g'] < 20) &
                (data_a['err_g'] < 0.1))
        data_a = data_a[mask]
        mask = shifts['field'] == f + '_1'
        s = shifts[mask]
        data_b = data_a.copy()
        for b in ['u', 'g', 'r1', 'r2', 'i', 'Ha']:
            data_b[b] = data_b[b] - s[b]

        tables_a.append(data_a)
        tables_b.append(data_b)

    t_a = vstack(tables_a)
    t_b = vstack(tables_b)

    t_b.write('outdata/carina/tables/carina_before_calib.fits')
    t_a.write('outdata/carina/tables/carina_after_calib.fits')


def cc_plots_all():

    t_a = Table.read('outdata/carina/tables/carina_after_calib.fits')
    t_b = Table.read('outdata/carina/tables/carina_before_calib.fits')

    cx = [('g', 'r2'), ('r1', 'i')]
    cy = [('u', 'g'), ('g', 'r2')]
    xlims = [[-0.5, 3], [-0.5, 2.5]]
    ylims = [[4.0, -2.0], [-0.5, 3.5]]
    mpl.rcParams['axes.linewidth'] = 0.1

    for i in range(len(cx)):
        plt.figure(figsize=(6, 2.5))
        ax1 = plt.subplot(1, 2, 1)
        binsx = np.arange(min(xlims[i]), max(xlims[i]) + 0.01, 0.01)
        binsy = np.arange(min(ylims[i]), max(ylims[i]) + 0.01, 0.01)
        ax2, vals = density_plot(t_b, cx[i], cy[i],
                                 0,
                                 0,
                                 ax1, [binsy, binsx])
        ax1.set_xlabel('${0} - {1}$'.format(cx[i][0][0], cx[i][1][0]),
                       fontsize=6)
        ax1.set_ylabel('${0} - {1}$'.format(cy[i][0][0], cy[i][1][0]),
                       fontsize=6)
        ax1.set_xlim(xlims[i])
        ax1.set_ylim(ylims[i])
        ax1.set_xticks(np.arange(xlims[i][0] + 0.5, xlims[i][1], 0.5))
        ax1.set_yticklabels(str(float(i)) for i in ax1.get_yticks())
        ax1.tick_params(labelsize=6, width=0.2, which="both")
        ax1.minorticks_on()
        ax2 = plt.subplot(1, 2, 2, sharey=ax1)
        ax2, vals = density_plot(t_a, cx[i], cy[i],
                                 0,
                                 0,
                                 ax2, [binsy, binsx])
        ax2.set_xlabel('${0} - {1}$'.format(cx[i][0][0], cx[i][1][0]),
                       fontsize=6)
        ax2.set_ylabel('${0} - {1}$'.format(cy[i][0][0], cy[i][1][0]),
                       fontsize=6)
        ax2.set_xticks(np.arange(xlims[i][0] + 0.5, xlims[i][1], 0.5))
        ax2.tick_params(labelsize=6, width=0.2, which="both")
        ax2.set_ylabel('')
        ax2.set_xlim(xlims[i])
        ax2.set_ylim(ylims[i])
        ax2.minorticks_on()
        plt.setp(ax2.get_yticklabels(), visible=False)
    # plt.tight_layout()
        plt.subplots_adjust(wspace=0.01)
        plt.savefig('outdata/calib/plots/all_{0}.pdf'.format(i),
                    bbox_inches='tight', rasterized=True)
###############################################################################


def merge(infile, outfile, conf_max=95):
    """

    :param infile:
    :param outfile:
    :param conf_max:
    """
    data = Table.read(infile)

    x, y, z = data['Av_conf_Ha_1'], data['Av_conf_Ha_2'], data['Av_conf_Ha_3']

    x_ha, y_ha, z_ha = data['Ha_1'], data['Ha_2'], data['Ha_3']

    x_err_ha, y_err_ha, z_err_ha = (data['err_Ha_1'],
                                    data['err_Ha_2'], data['err_Ha_3'])

    x = np.where(x != x, 0, x)

    y = np.where(y != y, 0, y)

    z = np.where(z != z, 0, z)

    new_conf = np.zeros(len(x))

    new_ha = np.zeros(len(x))

    new_err_ha = np.zeros(len(x))

    f1 = np.where((x >= conf_max) & (y >= conf_max) & (z >= conf_max))

    new_conf[f1] = (x[f1] + y[f1] + z[f1]) / 3.

    new_ha[f1] = (x_ha[f1] + y_ha[f1] + z_ha[f1]) / 3.

    new_err_ha[f1] = (x_err_ha[f1] + y_err_ha[f1] + z_err_ha[f1]) / 3.

    f2 = np.where((x >= conf_max) & (y >= conf_max) & (z <= conf_max))

    new_conf[f2] = (x[f2] + y[f2]) / 2.

    new_ha[f2] = (x_ha[f2] + y_ha[f2]) / 2.

    new_err_ha[f2] = (x_err_ha[f2] + y_err_ha[f2]) / 2.

    f3 = np.where((x <= conf_max) & (y >= conf_max) & (z >= conf_max))

    new_conf[f3] = (z[f3] + y[f3]) / 2.

    new_ha[f3] = (z_ha[f3] + y_ha[f3]) / 2.

    new_err_ha[f3] = (z_err_ha[f3] + y_err_ha[f3]) / 2.

    f4 = np.where((x >= conf_max) & (y <= conf_max) & (z >= conf_max))

    new_conf[f4] = (z[f4] + x[f4]) / 2.

    new_ha[f4] = (z_ha[f4] + x_ha[f4]) / 2.

    new_err_ha[f4] = (z_err_ha[f4] + x_err_ha[f4]) / 2.

    f5 = np.where((x >= conf_max) & (y <= conf_max) & (z <= conf_max))

    new_conf[f5] = x[f5]

    new_ha[f5] = x_ha[f5]

    new_err_ha[f5] = x_err_ha[f5]

    f6 = np.where((x <= conf_max) & (y >= conf_max) & (z <= conf_max))

    new_conf[f6] = y[f6]

    new_ha[f6] = y_ha[f6]

    new_err_ha[f6] = y_err_ha[f6]

    f7 = np.where((x <= conf_max) & (y <= conf_max) & (z >= conf_max))

    new_conf[f7] = z[f7]

    new_ha[f7] = z_ha[f7]

    new_err_ha[f7] = z_err_ha[f7]

    f8 = np.where((x <= conf_max) & (y <= conf_max) & (z <= conf_max))

    a = np.array([list(x[f8]), list(y[f8]), list(z[f8])])

    b = np.array([list(x_ha[f8]), list(y_ha[f8]), list(z_ha[f8])])

    c = np.array([list(x_err_ha[f8]), list(y_err_ha[f8]), list(z_err_ha[f8])])

    m = np.argmax(a, axis=0)

    new_conf[f8] = np.max(a, axis=0)

    new_ha[f8] = np.array([b[:, i][j] for i, j in enumerate(m)])

    new_err_ha[f8] = np.array([c[:, i][j] for i, j in enumerate(m)])

    data['Ha'] = new_ha
    data['err_Ha'] = new_err_ha
    data['Ha_conf'] = new_conf

    data.write(outfile, format='fits', overwrite=True)

###############################################################################


def sort(infile, outfile):

    data = Table.read(infile)

    bands = ['r1', 'r2', 'g', 'i', 'u']

    for b in bands:
        newm = np.zeros(len(data))
        newe = np.zeros(len(data))
        cols = fnmatch.filter(data.colnames, b + '_*')
        ecols = fnmatch.filter(data.colnames, 'err_' + b + '_*')

        delta = data[cols[0]] - data[cols[1]]
        mask1 = abs(delta) > 0.1
        mask2 = delta < 0
        mask3 = delta > 0

        mask4 = ((data[cols[0]] == data[cols[0]]) &
                 (data[cols[1]] != data[cols[1]]))
        mask5 = ((data[cols[0]] != data[cols[0]]) &
                 (data[cols[1]] == data[cols[1]]))

        mask6 = abs(delta) < 0.1

        mask7 = ((data[cols[0]] != data[cols[0]]) &
                 (data[cols[1]] != data[cols[1]]))

        newm[mask1 & mask2] = data[cols[0]][mask1 & mask2]
        newm[mask1 & mask3] = data[cols[1]][mask1 & mask3]
        newm[mask4] = data[cols[0]][mask4]
        newm[mask5] = data[cols[1]][mask5]

        newm[mask6] = np.mean([data[cols[0]][mask6],
                               data[cols[1]][mask6]], axis=0)

        newm[mask7] = np.repeat(np.nan, len(data[mask7]))

        newe[mask1 & mask2] = data[ecols[0]][mask1 & mask2]
        newe[mask1 & mask3] = data[ecols[1]][mask1 & mask3]
        newe[mask4] = data[ecols[0]][mask4]
        newe[mask5] = data[ecols[1]][mask5]
        newe[mask6] = np.sqrt([data[ecols[0]][mask6] ** 2,
                               data[ecols[1]][mask6] ** 2]) / 2.

        newe[mask7] = np.repeat(np.nan, len(data[mask7]))

        data[b] = newm
        data['err_' + b] = newe

    data.write(outfile, overwrite=True)

###############################################################################


def makeshifttable(median=False, concat=False):
    print '\nCompiling Shift Table...\n'
    files = glob.glob(wd + 'tables/*shifts.npy')
    files2 = glob.glob(wd + 'tables/*sigmas.npy')
    concats = []
    offsets = []
    u = []
    g = []
    r1 = []
    r2 = []
    i = []
    Ha = []
    u_s = []
    g_s = []
    r1_s = []
    r2_s = []
    i_s = []

    # if median:
    #     for fi in files:
    #         field = fi.split('vphas_')[1].split('-')[0]
    #         shifts = np.load(fi)[0]
    #         for o in [0, 1]:
    #             offsets.append(field + '_' + str(o + 1))
    #             u.append(np.mean([shifts['u'][0], shifts['u'][1]]))
    #             g.append(np.mean([shifts['g'][0], shifts['g'][1]]))
    #             r1.append(np.mean([shifts['r1'][0], shifts['r1'][1]]))
    #             r2.append(np.mean([shifts['r2'][0], shifts['r2'][1]]))
    #             i.append(np.mean([shifts['i'][0], shifts['i'][1]]))
    #             Ha.append(np.mean([shifts['r1'][0], shifts['r1'][1]]))

    concat_data = Table.read(wd + 'bin/vphas_infoblock_p97.fits')
    ld = []
    b = []
    groups = []
    offsets_t = []
    u_t = []
    g_t = []
    r1_t = []
    r2_t = []
    i_t = []
    Ha_t = []
    for fi in files:
        shifts = np.load(fi)[0]
        field = fi.split('vphas_')[1].split('-')[0]
        mask = concat_data['Field'] == 'vphas_' + field
        for o in [0, 1]:
            ld.append(concat_data['GAL_LONG'][mask][0])
            b.append(concat_data['GAL_LAT'][mask][0])
            groups.append(concat_data['GroupID'][mask][0])
            offsets_t.append(field + '_' + str(o + 1))
            u_t.append(shifts['u'][o])
            g_t.append(shifts['g'][o])
            r1_t.append(shifts['r1'][o])
            r2_t.append(shifts['r2'][o])
            i_t.append(shifts['i'][o])
            Ha_t.append(shifts['r1'][o])

    u_t = np.array(u_t)
    g_t = np.array(g_t)
    r1_t = np.array(r1_t)
    r2_t = np.array(r2_t)
    i_t = np.array(i_t)
    Ha_t = np.array(Ha_t)
    nexp = np.zeros(len(u_t))

    m = [u_t, g_t, r1_t, r2_t, i_t, Ha_t]

    for gr in np.unique(groups):
        mask = (groups == gr)
        for j, l in enumerate(m):
            if j == 0:
                mask2 = abs(m[1])[mask] < 0.2
            else:
                mask2 = abs(l)[mask] < 0.2
            if len(mask2[mask2]) == 0:
                l[mask] = np.median(l[mask])
            elif len(mask2[mask2]) == len(mask[mask]):
                l[mask] = np.median(l[mask])
            elif len(mask2[mask2]) > len(mask[mask]) / 2:
                l[mask] = np.median(l[mask][mask2])
            else:
                l[mask][mask2] = np.median(l[mask][mask2])
                l[mask][-mask2] = np.median(l[mask][-mask2])

    for fi in files:
        mask = concat_data['Field'] == fi.split('/')[-1].split('-')[0]
        field = fi.split('vphas_')[1].split('-')[0]
        shifts = np.load(fi)[0]
        for o in [0, 1]:
            # concats.append(concat_data['GroupID'][mask][0])
            # offsets.append(field + '_' + str(o + 1))
            u.append(shifts['u'][o])
            g.append(shifts['g'][o])
            r1.append(shifts['r1'][o])
            r2.append(shifts['r2'][o])
            i.append(shifts['i'][o])
            Ha.append(shifts['r1'][o])

    for fi in files2:
        shifts = np.load(fi)[0]
        for o in [0, 1]:
            # concats.append(concat_data['GroupID'][mask][0])
            # offsets.append(field + '_' + str(o + 1))
            #u_s.append(shifts['u'][o])
            g_s.append(shifts['g'][o])
            r1_s.append(shifts['r1'][o])
            r2_s.append(shifts['r2'][o])
            i_s.append(shifts['i'][o])
    print len(u_t), len(u)
    shift_table = Table([offsets_t, ld, b, u, g, r1, r2, i, Ha, u_t, g_t,
                         r1_t, r2_t, i_t, Ha_t, groups, nexp,
                         g_s, r1_s, r2_s, i_s],
                        names=['field', 'GAL_LONG', 'GAL_LAT', 'u', 'g', 'r1',
                               'r2', 'i', 'Ha', 'u_c', 'g_c', 'r1_c',
                               'r2_c', 'i_c', 'Ha_c', 'Concat', 'nexp',
                               'g_s', 'r1_s', 'r2_s', 'i_2'])
    shift_table.write(wd + 'tables/shifts_all.fits', overwrite=True)


def spatial_map_field(field):

    a2v = {'r1': 'rmag', 'r2': 'rmag', 'g': 'gmag', 'i': 'imag'}
    ab_to_vega = {'r1': -0.136, 'r2': -0.136, 'g': 0.123, 'i': -0.373}
    lims = {'r1': (13.5, 16.), 'r2': (13.5, 16.),
            'g': (14, 16.), 'i': (13.5, 15.),
            'u': (13, 20)}
    files = glob.glob('tables/vphas_{0}-apass.fits'.format(field))

    for f in files:
        data = Table.read(f)
        for offset in ['1', '2']:
            it = 0
            plt.figure(figsize=(12, 10))
            for band in ['r1', 'r2', 'g', 'i']:
                col = band + '_' + offset
                er = 'err_' + col
                mask1 = ((data[col] == data[col]) &
                         (data[a2v[band]] == data[a2v[band]]))

                mask2 = ((data[col] > lims[band][0]) &
                         (data[col] < lims[band][1]) &
                         (data[er] < 0.1) & (data[a2v[band]] < 99))
                mask = mask1 & mask2
                d = data[mask]
                delta = (d[band + '_' + offset] -
                         (d[a2v[band]] + ab_to_vega[band]))

                mask_delta = (np.abs(delta - np.median(delta)) >
                              (2 * np.std(delta)))
                ax = plt.subplot(221 + it)
                it += 1
                if len(mask_delta[mask_delta]) != 0:
                    a = ax.scatter(d['RA_' + offset][-mask_delta],
                                   d['Dec_' + offset][-mask_delta],
                                   c=delta[-mask_delta],
                                   edgecolor='None', s=25, cmap='jet')
                else:
                    a = ax.scatter(d['RA_' + offset],
                                   d['Dec_' + offset], c=delta,
                                   edgecolor='None', s=25, cmap='jet')
                ax.set_xlabel('RA')
                ax.set_ylabel('DEC')
                cb = plt.colorbar(a)
                cb.set_label(r'$\delta {0}$'.format(band), size=18)
                cb.solids.set_edgecolor("face")

            plt.tight_layout()
            plt.subplots_adjust(top=0.92, bottom=0.05)
            print 'Saving...'
            plt.suptitle(r'\underline{APASS Comparison Spatial Map}',
                         fontsize=20)
            plt.savefig('plots/spatial_map_{0}_{1}.png'.format(field, offset),
                        bbox_inches='tight', rasterized=True)
            plt.close()


def spatial_map(band, vmin=-0.25, vmax=0.25):

    a2v = {'r1': 'rmag', 'r2': 'rmag', 'g': 'gmag', 'i': 'imag'}
    ab_to_vega = {'r1': -0.136, 'r2': -0.136, 'g': 0.123, 'i': -0.373}
    lims = {'r1': (13.5, 16.), 'r2': (13.5, 16.),
            'g': (14, 16.), 'i': (13.5, 15.),
            'u': (13, 20)}
    offset = '2'

    data = Table.read('vphas_apass_all.fits')
    mask = (data['l'] > 0) & (data['l'] < 60)
    data['l'][mask] = data['l'][mask] + 360.
    col = band + '_' + offset
    er = 'err_' + col

    mask1 = (data[col] == data[col]) & (data[a2v[band]] == data[a2v[band]])

    mask2 = ((data[col] > lims[band][0]) & (data[col] < lims[band][1]) &
             (data[er] < 0.1) & (data[a2v[band]] < 99))
    mask = mask1 & mask2
    d = data[mask]
    delta = d[band + '_' + offset] - (d[a2v[band]] + ab_to_vega[band])
    plt.figure(figsize=(20, 20))
    lranges = [(400, 352.5), (352.5, 300.5), (300.5, 253,), (253, 205.5)]
    for it, lr in enumerate(lranges):
        mask = (d['l'] < lr[0]) & (d['l'] > lr[1])

        ax = plt.subplot(411 + it)
        a = ax.scatter(d['l'][mask], d['b'][mask], c=delta[mask],
                       edgecolor='None', s=2, vmin=vmin, vmax=vmax,
                       cmap='jet')
        ax.set_xlim(lr[0], lr[1])

        ax.set_xticks(np.arange(lr[0], lr[1], -20))
        ax.set_xticks(np.arange(lr[0], lr[1], -5), minor=True)
        labels = []
        for l in np.arange(lr[0], lr[1], -20):
            if l >= 360:
                labels.append(str(l - 360) + r'$^{\circ}$')
            else:
                labels.append(str(l) + r'$^{\circ}$')

        ax.set_xticklabels(labels, fontsize=24)
        ax.set_ylim(-5.9, 5.9)
        ax.set_yticklabels([str(int(i)) +
                            r'$^{\circ}$' for i in ax.get_yticks()],
                           fontsize=24)
        ax.set_ylabel('$b$', fontsize=24)
        ax.set_xlabel('$\ell$', fontsize=24)
        #ax.minorticks_on()
        ax.grid(linestyle='dashed', which='both', alpha=0.5)
        cb = plt.colorbar(a, pad=0.01)
        cb.set_label(r'$\delta {0}$'.format(band), size=24)
        cb.ax.tick_params(labelsize=20)
    print 'Saving..'
    plt.savefig('spatial_map_{0}_{1}_alt.png'.format(band, offset),
                bbox_inches='tight', rasterized=True, dpi=200)
    plt.close()


def table_stack(infits, outfits):
    """
    infits = list of input tables to be stacked
    outfits = name of output file

    :param infits:
    :param outfits:
    """

    command = stilts + " tcatn nin=" + str(len(infits))

    for i in range(len(infits)):
        command += " in" + str(i + 1) + "=" + infits[i]

    command += " out=" + outfits

    os.system(command)

def stack_apass_xmatches():
    files = glob.glob('tables/vphas_????-apass.fits')
    table_stack(files, 'vphas_apass_all.fits')


###############################################################################


if __name__ == "__main__":
    import sys
    os.chdir('/car-data/msmith/tools/calib/')
    args = sys.argv
    arg = args[1]
    print arg

    if args[-1] != 'apply' and args[-1] != 'makeshift' and args[-1] != 'spatial':
        if not os.path.exists(wd + 'tables/vphas_' + arg +
                              '-apass-shifts.npy'):
            f = vphas_field(arg)
            print '\nDownloading files...\n'
            ab = f.download()
            if ab != 'Missing files!':
                f.apass_calibrate()
            else:
                log = open('bin/{0}_log.txt'.format(arg), 'w')
                log.close()
        else:
            print 'Already Done!'

    elif args[-1] == 'apply':
        if not os.path.exists(wd + 'fields/vphas_' + arg + '.fits'):
            f = vphas_field(arg)
            print 'Applying Shifts...\n'
            f.aply_shifts()
            print 'Merging offsets...\n'
            f.sort()

    elif args[-1] == 'spatial':
        spatial_map_field(arg)
    else:
        print 'Already Done!'

    if args[-1] == 'makeshift':
        makeshifttable()
