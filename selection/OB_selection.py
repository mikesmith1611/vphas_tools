import sys
from astropy.io import fits
from astropy.table import vstack
import os
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.path as pltpath
import matplotlib.pyplot as plt
import vphas
glob = vphas.glob
np = vphas.np
fnmatch = vphas.fnmatch
Table = vphas.Table


if len(sys.argv) != 1:

    f = sys.argv[1]

else:
    f = raw_input('Which field? ')


def init(field):

    global in_field
    global out_plot
    global out_fits
    global root
    global stilts

    # Directories

    if os.uname()[1] == 'uhppc39.herts.ac.uk':
        in_field = 'outdata/calib/fields/vphas_' + field + '.fits'
        out_plot = 'outdata/selection/plots/OB_selection_' + field + '.pdf'
        out_fits = 'outdata/selection/tables/OB_selection_' + field + '.fits'
        root = '/local/home/msmith/vphas/'

    else:
        d = '/car-data/msmith/tools/'
        in_field = d + 'calib/fields/vphas_' + field + '.fits'
        out_plot = d + 'selection/plots/OB_selection_' + field + '.pdf'
        out_fits = d + 'selection/tables/OB_selection_' + field +\
                       '.fits'
        root = '/car-data/msmith/'
        stilts = 'java -jar /home/msmith/bin/topcat/topcat-*.jar -stilts'


in_known = 'None'

# Photometry Constraints

mag_limits = {'g': (13., 20.)}
err_limits = {'err_g': 0.1, 'err_r2': 0.1, 'err_i': 0.1}
class_exclusions = {'g': [0], 'r2': [0], 'i': [0]}
conf_cuts = {'u': 90, 'g': 90, 'r2': 90, 'i': 90}

# Reddening law used for selection

rv_low = [3.1, 3.8, 3.8, 3.8, 3.8]
rv_up = [3.8, 3.8, 3.8, 3.8, 3.8]

# Spectral types to make selection

spectypes = ['B3_V', 'B1_V', 'R--J']
u_g_shift = [0.1, 0, 0, 0, 0]

# Plot options

spectype_colors = ['b'] * 3
spectype_markers = ['+'] * 3

redline_colours = ['k'] * 3

marker_labels = ['B3V-B2V', 'B2V-B1V', 'B1V-O6V', 'O6V-Earlier', 'R-J-Above']

known_object_labels = ['Known Members']

###############################################################################


def OBplot(OBstars, background, spectypes, rv_low, edgecolors,
           redline_colours, markers,
           labels,
           filename, u_g_shift, known_objects='None'):
    """

    :param OBcut:
    :param background:
    :param spectypes:
    :param rv_low:
    :param rv_up:
    :param edgecolors:
    :param markers:
    :param labels:
    :param filename:
    :param u_g_shift:
    :param known_objects:
    """
    print 'Plotting C-C diagram ....\n'

    fig = plt.figure(figsize=(6, 6))
    ax1 = fig.add_axes([0., 0., 1., 1.])

    ua, ga, ra = [background['u'], background['g'], background['r2']]

    ax1.plot(ga - ra, ua - ga, marker='.',
             linestyle='None', linewidth=0, mec='#888888',
             mfc='#888888', ms=1, alpha=0.2, rasterized=True)

    for i in range(len(OBstars)):
        u, g, r = [OBstars[i]['u'], OBstars[i]['g'], OBstars[i]['r2']]
        ax1.plot(g - r, u - g,
                 linestyle='None', marker=markers[i], mfc=edgecolors[i],
                 mec=edgecolors[i],
                 label=labels[i], ms=5, alpha=0.5, rasterized=False)

    ZAMS = getmodelcolours(ms=True)

    mask = np.argsort(ZAMS['g'] - ZAMS['r'])

    ZAMS = ZAMS[mask]

    ax1.plot(ZAMS['g'] - ZAMS['r'],
             ZAMS['u'] - ZAMS['g'],
             marker='None', linestyle='-', color='k', ms=2, rasterized=False)

    for i in range(len(spectypes)):
        model = getmodelcolours(spectype=spectypes[i], R_V=rv_low[i])

        ax1.plot(model['g'][0::10] - model['r'][0::10],
                 model['u'][0::10] - model['g'][0::10] + u_g_shift[
                     i],
                 marker='x', linestyle='-', color=redline_colours[i],
                 ms=5, rasterized=False)

    plt.xlabel(r'$g-r$', fontsize=20)
    plt.ylabel(r'$u-g$', fontsize=20)
    plt.xlim([-0.6, 2.75])
    plt.ylim([-1.8, 2.0])
    # h,l = ax1.get_legend_handles_labels()
    # mask = np.array([4,3,2,1,0,5],dtype=int)
    # h,l = np.array(h)[mask],np.array(l)[mask]
    # legend = ax1.legend(h,l,scatterpoints = 1,
    #   numpoints=1,fontsize=14,loc='lower left',markerscale=2)
    # frame = legend.get_frame()
    # frame.set_linewidth(0.5)
    ax1.text(-0.11 - 0.1, 0, 'MS', rotation=0, fontsize=14)
    ax1.text(-0.56 + 0.07 - 0.1, -1.67, 'R-J', rotation=0, fontsize=14,
             color=redline_colours[2])
    # ax1.text(-0.54 + 0.07 - 0.1, -1.48, 'O6V', rotation=0, fontsize=14,
    #        color=edgecolors[1])
    ax1.text(-0.48 + 0.08 - 0.1, -1.33, 'B1V', rotation=0, fontsize=14,
             color=redline_colours[2])
    # ax1.text(-0.46 + 0.08 - 0.1, -1.22, 'B2V', rotation=0, fontsize=14,
    #         color=edgecolors[1])
    ax1.text(-0.44 + 0.09 - 0.1, -1.05, 'B3V', rotation=0, fontsize=14,
             color=redline_colours[0])
    ax1.text(0.08, -0.53, '$1A_0$', rotation=0, fontsize=14)

    plt.gca().invert_yaxis()

    for s in ['G0_V']:
        model_GOV = getmodelcolours(spectype=s)
        GOV_gr = model_GOV['g'] - model_GOV['r']
        GOV_ug = model_GOV['u'] - model_GOV['g']
        # GOV_ri = model_GOV.field('r') - model_GOV.field('i')
        # GOV_rha = model_GOV.field('r') - model_GOV.field('ha')
        ax1.plot(GOV_gr[0::10], GOV_ug[0::10],
                 marker='x', linestyle='-', color='k', ms=5)
    # ax1.text(1.76-0.1,1.2,'$R_V = 3.1$',rotation=-46,fontsize=14)
    ax1.text(1.2, 0.8, 'G0V', rotation=0, fontsize=14)

    print 'Saving figure ....'
    plt.savefig(filename, bbox_inches='tight')
    plt.close()

###############################################################################


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


def OBcut(data, spectypes, **kwargs):
    """
    ====================================================================
    OBcut - makes a colour cut above the specified reddening vector
    ====================================================================
    Variables

    redvectorfilename - file path of reddening vector

    data - fits record array of data

    **kwargs

    shift - u band shift

    :param data:
    :param spectypes:
    :param kwargs:
    """
    d = {'ushift': 0, 'gshift': 0, 'R_V_low': [3.1], 'R_V_up': [5.0],
         'u_g_shift': [0]}
    for kwarg in kwargs:
        d[kwarg] = kwargs[kwarg]

    out = []

    for i in range(len(spectypes)):

        reddata1 = getmodelcolours(R_V=str(d['R_V_low'][i]),
                                   spectype=spectypes[i])

        z = np.polyfit(reddata1['g'] - reddata1['r'],
                       reddata1['u'] - reddata1['g'] + d['u_g_shift'][i], 3)

        mask2 = ((data['u'] - data['g']) <
                 np.lib.polyval(z, data['g'] - data['r2']))

        if i + 1 < len(spectypes):
            reddata2 = getmodelcolours(R_V=str(d['R_V_up'][i]),
                                       spectype=spectypes[i + 1])
            z2 = np.polyfit(reddata2['g'] - reddata2['r'],
                            reddata2['u'] - reddata2['g'] +
                            d['u_g_shift'][i + 1], 3)
            mask4 = ((data['u'] - data['g']) >
                     np.lib.polyval(z2, data['g'] - data['r2']))

            cut = data[mask2 & mask4]

            out.append(cut)
        else:
            cut = data[mask2]

            out.append(cut)
    return out


###############################################################################

def data_constraints(data, mag_limits=None, err_limits=None,
                     class_exclusions=None, conf_cuts=None):
    """

    :param data:
    :param mag_limits:
    :param err_limits:
    :param class_exclusions:
    :return:
    """
    if mag_limits is not None:
        for key in mag_limits.keys():
            mask = ((data[key] >= mag_limits[key][0]) &
                    (data[key] <= mag_limits[key][1]))
            data = data[mask]

    if err_limits is not None:
        for key in err_limits.keys():
            mask = (data[key] <= err_limits[key])
            data = data[mask]

    if class_exclusions is not None:
        for key in class_exclusions.keys():
            x = class_exclusions[key]
            mask = ((data['Class_' + key + '_1'] != x) |
                    (data['Class_' + key + '_2'] != x))

            data = data[mask]

    if conf_cuts is not None:
        for key in conf_cuts.keys():

            mask4 = ((data[key + '_1'] == data[key + '_1']) &
                     (data[key + '_2'] != data[key + '_2']) &
                     (data['Av_conf_' + key + '_1'] < conf_cuts[key]))

            mask5 = ((data[key + '_1'] != data[key + '_1']) &
                     (data[key + '_2'] == data[key + '_2']) &
                     (data['Av_conf_' + key + '_2'] < conf_cuts[key]))

            data = data[-mask4 & -mask5]

    return data


###############################################################################


def OB_selction():
    """

    :param OB_selection_cfg:
    """

    print 'Loading data ....\n'
    print in_field

    data = Table.read(in_field)

    # mask = (data['logteff_median'] > 4.3)

    # data = data[mask]

    print 'Applying constraints....\n'

    data = data_constraints(data,
                            mag_limits=mag_limits,
                            err_limits=err_limits,
                            class_exclusions=class_exclusions,
                            conf_cuts=conf_cuts)

    print 'Making OB cut....\n'

    OBcuts = OBcut(data, spectypes,
                   ushift=0,
                   gshift=0,
                   R_V_low=rv_low,
                   R_V_up=rv_up,
                   u_g_shift=u_g_shift)

    print 'Saving data files....\n'

    vstack(OBcuts).write(out_fits, overwrite=True)

    OBplot(OBcuts, data, spectypes,
           rv_low,
           spectype_colors,
           redline_colours,
           spectype_markers,
           marker_labels,
           out_plot, u_g_shift)
    return True
    # except:
    #    print 'No data for ' + f
    #    return False


###############################################################################


def xmatch_2mass(infile, outfile, nir_cat='2MASS'):
    """

    :param infile:
    :param outfile:
    :param nir_cat:
    """

    cats = {'2MASS': 'II/246/out', 'Ascenso': 'J/A+A/466/137/w2phot'}

    os.system(stilts + ' cdsskymatch cdstable=' + cats[nir_cat] +
              ' find=best ' +
              'in=' + infile + ' ra=RA dec=DEC radius=1 out=' + outfile)


###############################################################################

def all_params(infile):
    """

    :param infile:
    """
    data = Table.read(infile)

    mask = ((data['r2'] > 0) & (data['g'] > 0) &
            (data['i'] > 0) & (data['u'] > 0) &
            (data['err_r2'] > 0) & (data['err_g'] > 0) &
            (data['err_i'] > 0) & (data['err_u'] > 0) &
            (data['Jmag'] > 0) & (data['Hmag'] > 0) &
            (data['Kmag'] > 0) & (data['e_Jmag'] > 0) &
            (data['e_Hmag'] > 0) & (data['e_Kmag'] > 0))

    data = data[mask]

    data.write(infile, overwrite=True)

###############################################################################


def make_sed_pbs():
    nodes = '20'
    template = open('../sed/mcmc.pbs', 'r+')
    newfname = '../sed/pbs/' + f + '_mcmc.pbs'
    os.system('rm ' + newfname)
    newfile = open(newfname, 'w')

    for line in template:

        if line.startswith('#PBS -N'):
            newfile.write('#PBS -N sedmcmc_' + f + '\n')

        elif line.startswith(r'#PBS -l nodes=100:ppn=1'):
            newfile.write(r'#PBS -l nodes=' + nodes + ':ppn=1\n')

        elif line.startswith(r'#PBS -l pmem=1gb,walltime=2:00:00'):

            t = (len(Table.read(out_fits[:-5] + '_2mass.fits')) /
                 float(nodes)) * 60

            t1 = int(60 * np.ceil(t * 2 / 60))

            m, s = divmod(t1, 60)
            h, m = divmod(m, 50)
            wt = "%d:%02d:%02d" % (h, m, s)

            print wt
            newfile.write(r'#PBS -l pmem=1gb,walltime=' + wt + '\n')

        elif line.startswith('python'):
            newfile.write('python /car-data/msmith/tools/sed/sed_mcmc.py ' +
                          f + '\n')

        elif line.startswith('    ipcluster'):
            newfile.write(line.replace('sedmcmc', 'sedmcmc_' + f))

        elif line.startswith('ipcluster'):
            newfile.write(line.replace('sedmcmc', 'sedmcmc_' + f))

        elif line.startswith('while'):
            newfile.write(line.replace('sedmcmc', 'sedmcmc_' + f))

        elif line.startswith('if'):
            newfile.write(line.replace('sedmcmc', 'sedmcmc_' + f))

        else:
            newfile.write(line)

    newfile.close()
    template.close()

###############################################################################


if __name__ == "__main__":
    import sys
    os.chdir('/car-data/msmith/tools/selection/')
    init(f)

    f = sys.argv[1]

    a = OB_selction()

    if a:
        print out_fits

        xmatch_2mass(out_fits, out_fits[:-5] + '_2mass.fits')

        all_params(out_fits[:-5] + '_2mass.fits')

        dat = Table.read(out_fits[:-5] + '_2mass.fits')

        dat['ID'] = np.arange(1, len(dat) + 1)

        dat['Field'] = np.array([f] * len(dat))

        dat.write(out_fits[:-5] + '_2mass.fits', overwrite=True)

        make_sed_pbs()
