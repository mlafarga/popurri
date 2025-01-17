"""
Build a telluric mask from a spectrum obtained with the ESO Sky Model Calculator (https://www.eso.org/observing/etc/skycalc Noll et al. 2012).

The installed script is called `tellmask`.

Example run:
tellmask PATH/TO/MODEL.fits --wmin 3600 --wmax 7100 --fk trans_ma --pltext pdf --pltsv --verbose
"""
import argparse
import os
import sys
import textwrap
import ipdb

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
# from scipy.signal import find_peaks, argrelmin, argrelmax

from popurri import peakutils
from popurri import plotutils
from popurri.plotutils import wavelength_label
from popurri import telluricutils

# mpl.rcdefaults()
plotutils.mpl_custom_basic()
plotutils.mpl_size_same(font_size=18)

# Strings
angstromstr = r'$\mathrm{\AA}$'
angstromstrunit = r'$[\mathrm{\AA}]$'

###############################################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description=textwrap.dedent('''
        `tellmask_from_esoskycalc.py`

        Build a telluric mask from a spectrum obtained with the ESO Sky Model Calculator (https://www.eso.org/observing/etc/skycalc Noll et al. 2012).
        '''),
        epilog=textwrap.dedent('''
            '''),
        )

    # Template
    parser.add_argument('filtpl', help='Input spectrum file.', type=str)
    parser.add_argument('--tpltype',  help='Input model type.', choices=['esoskycalc'],type=str, default='esoskycalc')
    parser.add_argument('--fk', help='Flux key to be used.', choices=['ftransmission', 'trans_o3', 'trans_ma', 'trans_rs', 'trans_ms'], default='ftransmission')

    # Mask parameters
    parser.add_argument('--wmin', help='Minimum wavelength of the mask [A]. If None, use minimum wavelength of the template provided.', type=float, default=None)
    parser.add_argument('--wmax', help='Maximum wavelength of the mask [A]. If None, use maximum wavelength of the template provided.', type=float, default=None)
    parser.add_argument('--lisfcut', help='Flux threshold(s) to select telluric lines (0 to 1). I.e. fgood >= fcut, or fbad < fcut.', nargs='+', type=float, default=[0.9, 0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999])
    parser.add_argument('--dw', help='Line "limit" to make a 0110 mask [A].', type=float, default=1e-8)

    # Initial cuts to ease visualisation
    parser.add_argument('--npixmin', help='', type=int, default=3)
    parser.add_argument('--depthmin', help='', type=float, default=1-5*1e-5)

    # Outputs
    parser.add_argument('--dirout', help='General output directory. Absolute path.', type=str, default='./tellmask_output/')
    parser.add_argument('--diroutplot', help='Directory inside DIROUT where to store general plots (DIROUT/DIROUTPLOTS/)', type=str, default='plots/')
    parser.add_argument('--filoutmask', help='Name of the mask file. If None (default), name of the input spectrum template changing extension to `.mas`.', type=str, default=None)

    # Plots
    parser.add_argument('--cmap', help='Colormap to be used for the masks plot.', type=str, default='plasma')
    parser.add_argument('--pltsv', help='Make and save plots.', action='store_true')
    parser.add_argument('--pltsh', help='Show plots.', action='store_true')
    parser.add_argument('--pltext', nargs='+', help='Extensions of the plots to be saved (e.g. `--pltext pdf png`)', default='pdf')

    # Script behaviour
    parser.add_argument('--verbose', help='', action='store_true')
    parser.add_argument('--pause', help='Pause inside python after running.', action='store_true')

    args = parser.parse_args()
    return args


def main():

    args = parse_args()

    verboseprint = print if args.verbose else lambda *a, **k: None

    verboseprint('\n')
    verboseprint('#' * 50)
    verboseprint('\nCreate telluric mask from a telluric model\n')
    verboseprint('#' * 50)
    verboseprint('\n')

    verboseprint(f'Input spectrum model: {args.filtpl} ({args.tpltype}, using {args.fk} flux key)')
    verboseprint(f'Include tellurics below: {args.lisfcut}')
    verboseprint(f'Mask w limits: {args.wmin} to {args.wmax} A')
    verboseprint('Masks saved in:', args.dirout)
    verboseprint('Plots saved in:', args.diroutplot)

    # Extend directories and files
    if isinstance(args.filtpl, str): args.filtpl = os.path.expanduser(args.filtpl)
    if isinstance(args.dirout, str): args.dirout = os.path.expanduser(args.dirout)
    if isinstance(args.diroutplot, str): args.diroutplot = os.path.expanduser(args.diroutplot)

    # Join directories
    args.diroutplot = os.path.join(args.dirout, args.diroutplot)

    # Create output directories if they do not exist
    for d in [args.dirout, args.diroutplot]:
        if not os.path.exists(d): os.makedirs(d)

    # To plot or not to plot
    doplot = args.pltsh or args.pltsv

    # Make sure figure extensions is a list
    if not isinstance(args.pltext, list): args.pltext = [args.pltext]

    # File ID
    filtplid = os.path.splitext(os.path.basename(args.filtpl))[0]

    # Flux cut values from less to more restrictive
    args.lisfcut = np.sort(args.lisfcut)[::-1]
    
    # Default wavelength keyword
    wk = 'w'

    # Wavelength unit
    xunit = 'A'

    ###########################################################################

    # Read input spectrum
    verboseprint('Reading input spectrum', args.filtpl)
    if args.tpltype == 'esoskycalc':
        model = telluricutils.ModelEsoSkyCalc(args.filtpl, wk=wk, fk=args.fk, tag=filtplid, dirout=args.dirout, diroutplot=args.diroutplot)
    else:
        sys.exit('File type not valid.')
    
    # Cut wavelength
    verboseprint(f'Cutting w from {args.wmin} to {args.wmax} A')
    model.cut_w(args.wmin, args.wmax)

    # Plot initial spectrum (after wavelength cut)
    if doplot:
        verboseprint('Plotting spectrum')
        # Plot spectrum (selected component fk)
        model.fig_spec(lisfk=[args.fk], sh=args.pltsh, sv=args.pltsv, svext=args.pltext)
        # Plot all transmission spectrum components
        model.fig_spec_transmission(sh=args.pltsh, sv=args.pltsv, svext=args.pltext)
    

    ###########################################################################

    # # Normalise ??????????????
    # if args.trans_ma is False:
    #     sys.exit('Normalisation not implemented yet.')

    ###########################################################################

    # Find lines
    verboseprint('Finding lines')
    model.find_lines()

    # Initial line cleaning
    # --- Remove lines with less than npixmin pixels
    verboseprint(f'Plot: Removing lines with less than {args.npixmin} pixels')
    model.rmv_lines_npix(npixmin=args.npixmin)
    # --- Remove very shallow lines (issues with the template)
    # ------ Only to reduce the number of lins to be plotted. Requirement with fcut will be stronger on depth anyway.
    verboseprint(f'Plot: Removing lines with more than {args.depthmin} depth')
    model.rmv_lines_depth(depthmin=args.depthmin)

    if doplot:
        verboseprint('Plotting initial lines')
        model.fig_spec_lines(sh=args.pltsh, sv=args.pltsv, svext=args.pltext)

    ###########################################################################

    # Make telluric mask
    verboseprint('Making mask including lines with more than {args.lisfcut} depth')
    for j, fcut in enumerate(args.lisfcut):  # for different fcut values
        model.make_mask(fcut=fcut, dw=args.dw, svobj=True)


    ###########################################################################

    # Plot the different masks
    verboseprint('Plotting masks')

    # Number of masks
    nmask = len(model.masks.keys())

    # Color per mask from colormap
    cmap = mpl.cm.get_cmap(args.cmap)
    norm = mpl.colors.Normalize(vmin=0, vmax=nmask - 1)

    # -------------------------------------------

    # One mask per panel
    fig, ax = plt.subplots(nmask, 1, figsize=(16, 2.5 * nmask), sharex=True, sharey=True, constrained_layout=True)
    for i, fcut in enumerate(args.lisfcut):
        # Plot spec
        ax[i] = model.plot_spec(ax=ax[i], xlabel='', ylabel='', rasterized=True)
        if i == 0: xlim, ylim = ax[i].get_xlim(), ax[i].get_ylim()
        # Plot mask
        ax[i] = model.masks[fcut].plot_mask(ax=ax[i], fscale=ylim[1],  xlabel='', ylabel='', zorder=0, leglabel=f'fcut={fcut:.6f}', color='0.5', alpha=0.6)
        ax[i].set(xlim=xlim, ylim=ylim)
    # Format
    for a in ax.flatten():
        a.minorticks_on()
        a.set_ylabel('Flux')
        a.legend(fontsize='small', loc='lower left')
    ax[0].set(xlim=xlim, ylim=ylim)
    ax[-1].set_xlabel(wavelength_label(x=xunit))
    ax[0].set_title(filtplid, loc='left')
    plotutils.figout(fig, filout=os.path.join(args.diroutplot, filtplid + '_masks_panels'), sv=args.pltsv, svext=args.pltext, sh=args.pltsh)

    # -------------------------------------------

    # One mask per figure
    for i, fcut in enumerate(args.lisfcut):
        fig, ax = plt.subplots(1, 1, figsize=(16, 4), constrained_layout=True)
        # Plot spec
        ax = model.plot_spec(ax=ax, xlabel='', ylabel='', rasterized=True)
        # Plot mask
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        ax = model.masks[fcut].plot_mask(ax=ax, fscale=ylim[1],  xlabel='', ylabel='', zorder=0, color='0.5', alpha=0.6)
        ax.set(xlim=xlim, ylim=ylim)
        # Format
        ax.minorticks_on()
        ax.set(ylabel='Flux')
        ax.set_xlabel(wavelength_label(x=xunit))
        ax.set_title(filtplid + '\n' + f'Mask fcut={fcut:.6f}', loc='left')
        plotutils.figout(fig, filout=os.path.join(args.diroutplot, filtplid + '_mask_fcut{:.6f}'.format(fcut)), sv=args.pltsv, svext=args.pltext, sh=args.pltsh)

    # -------------------------------------------

    # All masks same panel, different heights
    fig, ax = plt.subplots(1, 1, figsize=(16, 4), constrained_layout=True)
    # Plot spec
    ax = model.plot_spec(ax=ax, xlabel='', ylabel='', rasterized=True, zorder=nmask * 10)
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    for i, fcut in enumerate(args.lisfcut):
        ax = model.masks[fcut].plot_mask(ax=ax, fscale=ylim[1]-i/10,  xlabel='', ylabel='', zorder=0,  leglabel=f'fcut={fcut:.6f}', color=cmap(norm(i)), alpha=0.4)
    ax.minorticks_on()
    ax.set(xlim=xlim, ylim=ylim, xlabel=wavelength_label(x=xunit), ylabel='Flux')
    ax.legend(fontsize='small', loc='lower left')
    ax.set_title(filtplid, loc='left')
    plotutils.figout(fig, filout=os.path.join(args.diroutplot, filtplid + '_masks'), sv=args.pltsv, svext=args.pltext, sh=args.pltsh, cl=False)

    # Same as above but zoom flux to top 1%, below maximum of spec flux
    # Masks go in front of the spectrum here, to be able to see them (if not spectrum covers them in strong absorption regions)
    fig, ax = plt.subplots(1, 1, figsize=(16, 4), constrained_layout=True)
    yzoom = 0.01
    ylim1 = model.data[args.fk].max()
    dy = ylim1 - ylim[0]
    ylim0 = ylim1 - dy * yzoom
    # Plot spec
    ax = model.plot_spec(ax=ax, xlabel='', ylabel='', rasterized=True, zorder=-100)
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    for i, fcut in enumerate(args.lisfcut):
        ax = model.masks[fcut].plot_mask(ax=ax, fscale=ylim1-i/(10/yzoom),  xlabel='', ylabel='', zorder=0,  leglabel=f'fcut={fcut:.6f}', color=cmap(norm(i)), alpha=0.4)
    ax.minorticks_on()
    # ylim0 = ylim[0] + dy * (1 - yzoom)
    ax.set(xlim=xlim, ylim=(ylim0, ylim1), xlabel=wavelength_label(x=xunit), ylabel='Flux')
    ax.legend(fontsize='small', loc='lower left')
    ax.set_title(filtplid, loc='left')
    plotutils.figout(fig, filout=os.path.join(args.diroutplot, filtplid + '_masks_zoomy1'), sv=args.pltsv, svext=args.pltext, sh=args.pltsh, cl=False)

    ###########################################################################

    if args.pause: ipdb.set_trace()

    return

###############################################################################


if __name__ == "__main__":

    print('Running:', sys.argv)

    # main(sys.argv[1:])
    main()

    print('End')
