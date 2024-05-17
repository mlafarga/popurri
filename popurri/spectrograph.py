"""
Spectrograph-related functions

TODO CRIRES+ plotting different detectons
TODO ESPRESSO, MAROON-X plotting different slices
"""
import os
import sys

import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

from .plotutils import wavelength_label

dirhere = os.path.dirname(__file__)
dirdata = os.path.join(dirhere, './data/')

###############################################################################

# General plotting and strings

dictprop_nice = {
    'inst': 'Instrument str',
    'R': 'R',
    'spectral_sampling_px': 'Mean sampling',
    'pixel_ms': 'Pixel sampling [m/s]',
    'wmin_nm': '$\lambda_\mathrm{min}$ [nm]',
    'wmax_nm': '$\lambda_\mathrm{max}$ [nm]',
    'ndet': 'Number of detectors',  # in the reduced data
    'ndet_raw': 'Real number of detectors',
    'nord': 'Number of orders',
    'ord_ref': 'Reference order',
    'nslice': 'Number of slices',  # in the reduced data
    'nslice_raw': 'Real number of slices',
    'npix': 'Number of pixels spectral direction',
    'ron': 'Readout noise [e-]',
    'conversion_e_adu': 'Conversion factor [e-/ADU]',
    'position': 'Position',
    'type': 'Type',
    'inst_nice': 'Instrument',
    'inst_acronym_nice': 'Instrument acronym',
    'tel': 'Telescope str',
    'tel_nice': 'Telescope',
    'tel_acronym_nice': 'Telescope acronym',
    'tel_diameter': 'Telescope diameter',
    'observatory': 'Observatory str',
    'observatory_nice': 'Observatory',
    'observatory_acronym_nice': 'Observatory acronym',
    'year_start': 'Year start',
    'year_end': 'Year end',
    'ref': 'Reference',
    'notes': 'Notes',
}


###############################################################################

# Spectrograph

filprop = os.path.join(dirdata, 'spectrograph/spectrograph_properties.csv')


def read_spectrograph_properties(filprop=filprop):
    """
    Read spectrograph properties from file
    """
    df = pd.read_csv(filprop).set_index('inst')
    return df

###############################################################################

# Orders

def get_filords(inst):
    filords = os.path.join(dirdata, f'spectrograph/{inst}_orders.csv')
    return filords


def read_ords_carmvis(filords):
    ords_carmvis = pd.read_csv(filords, comment='#', skipinitialspace=True)
    # Transform from A to nm
    ords_carmvis['wmin_nm'] = ords_carmvis['wmin_A'] / 10.
    ords_carmvis['wmax_nm'] = ords_carmvis['wmax_A'] / 10.
    return ords_carmvis


def read_ords_carmnir(filords, ordcut=True):
    ords_carmnir = pd.read_csv(filords, comment='#', skipinitialspace=True)
    # Transform from A to nm
    ords_carmnir['wmin_nm'] = ords_carmnir['wmin_A'] / 10.
    ords_carmnir['wmax_nm'] = ords_carmnir['wmax_A'] / 10.
    # Cut orders
    if ordcut:
        'ordcut1'
        df_new = ords_carmnir.copy()
        df_new.set_index('ordcut1', inplace=True, drop=False)
        df_new2 = ords_carmnir.copy()
        df_new2.set_index('ordcut2', inplace=True, drop=False)
        pd.concat([df_new, df_new2]).sort_index()
        ords_carmnir = pd.concat([df_new, df_new2]).sort_index()
    return ords_carmnir


def read_ords_crires(filords):
    # TODO
    pass


def read_ords_criresplus(filords):
    ords = pd.read_csv(filords, comment='#', skipinitialspace=True).set_index(['setting', 'order'])

    # Transform from nm to A
    ords['wmin_det1_A'] = ords['wmin_det1_nm'] * 10.
    ords['wmax_det1_A'] = ords['wmax_det1_nm'] * 10.
    ords['wmin_det2_A'] = ords['wmin_det2_nm'] * 10.
    ords['wmax_det2_A'] = ords['wmax_det2_nm'] * 10.
    ords['wmin_det3_A'] = ords['wmin_det3_nm'] * 10.
    ords['wmax_det3_A'] = ords['wmax_det3_nm'] * 10.
    # Transform from nm to micron
    ords['wmin_det1_mu'] = ords['wmin_det1_nm'] * 1.e-3
    ords['wmax_det1_mu'] = ords['wmax_det1_nm'] * 1.e-3
    ords['wmin_det2_mu'] = ords['wmin_det2_nm'] * 1.e-3
    ords['wmax_det2_mu'] = ords['wmax_det2_nm'] * 1.e-3
    ords['wmin_det3_mu'] = ords['wmin_det3_nm'] * 1.e-3
    ords['wmax_det3_mu'] = ords['wmax_det3_nm'] * 1.e-3
    return ords


def get_ghosts_criresplus():
    inst = 'criresplus'
    ghosts = pd.read_csv(os.path.join(dirin, '{}_ghosts.csv'.format(inst)), comment='#', skipinitialspace=True)

    # Transform from nm to A
    ghosts['wmin_A'] = ghosts['wmin_nm'] * 10.
    ghosts['wmax_A'] = ghosts['wmax_nm'] * 10.
    # Transform from nm to micron
    ghosts['wmin_mu'] = ghosts['wmin_nm'] * 1.e-3
    ghosts['wmax_mu'] = ghosts['wmax_nm'] * 1.e-3
    return ghosts


# espresso_uhr11
# espresso_hr11
# espresso_hr21
# espresso_hr42
# espresso_mr42
# espresso_mr84

def read_ords_espresso(filords):
    # Assume HR mode
    return read_ords_espresso_hr(filords)


def read_ords_espresso_hr(filords):
    """
    """
    ords_espresso = pd.read_csv(filords, comment='#', skipinitialspace=True)
    # Transform from nm to A
    ords_espresso['wmin_A'] = ords_espresso['wmin_nm'] * 10.
    ords_espresso['wmax_A'] = ords_espresso['wmax_nm'] * 10.

    # Column with slices as strings
    ords_espresso['slices_str'] = ['{}, {}'.format(ords_espresso.loc[o]['ord_hr1'] - 1, ords_espresso.loc[o]['ord_hr2'] - 1) for o in ords_espresso.index]

    # Get slices instead of orders (HR, UHR)
    # i.e. duplicare rows but row_index refers to the slice, not the order
    df_new = ords_espresso.copy()
    df_new.set_index('ord_hr1', inplace=True, drop=False)
    df_new2 = ords_espresso.copy()
    df_new2.set_index('ord_hr2', inplace=True, drop=False)
    pd.concat([df_new, df_new2]).sort_index()
    ords_espresso_hr = pd.concat([df_new, df_new2]).sort_index()
    
    # Start order index at 0 instead of 1
    ords_espresso_hr['ord'] = ords_espresso_hr.index - 1
    ords_espresso_hr.set_index('ord', inplace=True, drop=False)

    return ords_espresso_hr


def read_ords_espresso_mr(filords):
    """
    """
    inst = 'espresso'
    ords_espresso = pd.read_csv(filords, comment='#', skipinitialspace=True)
    # Transform from nm to A
    ords_espresso['wmin_A'] = ords_espresso['wmin_nm'] * 10.
    ords_espresso['wmax_A'] = ords_espresso['wmax_nm'] * 10.

    # Start order index at 0 instead of 1, MR
    ords_espresso_mr = ords_espresso.copy()
    ords_espresso_mr['ord'] = ords_espresso_mr['ord_mr'] - 1
    ords_espresso_mr.set_index('ord', inplace=True, drop=False)
    ords_espresso_mr.sort_index(inplace=True)
    return ords_espresso_mr
    

def read_ords_expres(filords):
    # TODO
    pass


def read_ords_harps(filords):
    # inst = 'harps'
    # dataord = pd.read_csv(os.path.join(dirin, '{}_orders.csv'.format(inst)), comment='#', skipinitialspace=True, index_col=0)
    dataord = pd.read_csv(filords, comment='#', skipinitialspace=True, index_col=0)
    dataord.sort_index(inplace=True)
    dataord['ord'] = dataord.index

    # Transform from nm to A
    dataord['wmin_A'] = dataord['wmin_nm'] * 10.
    dataord['wmax_A'] = dataord['wmax_nm'] * 10.
    return dataord


def read_ords_harpsn(filords):
    # TODO
    pass


def read_ords_maroonx(filords):
    # TODO
    pass


def read_ords_neid(filords):
    # TODO
    pass


def read_ords(inst, filords, **kwargs):
    """
    Read orders of specific instrument.
    """
    # dataord = pd.read_csv(filords)
    # if inst == 'harps':
    #     dataord = read_ords_harps(filords)
    dictinst = {
        'carmvis': read_ords_carmvis,
        'carmnir': read_ords_carmnir,
        'crires': read_ords_crires,
        'criresplus': read_ords_criresplus,
        'espresso': read_ords_espresso,
        'espresso_hr': read_ords_espresso_hr,
        'express': read_ords_expres,
        'harps': read_ords_harps,
        'harpsn': read_ords_harpsn,
        'maroonx': read_ords_maroonx,
        'neid': read_ords_neid,
    }
    # Remove extra kwargs
    if inst != 'carmnir':
        kwargs.pop('ordcut')
    try:
        dataord = dictinst[inst](filords, **kwargs)
    except:
        print(f'Instrument {inst} not implemented')
    return dataord


###############################################################################

class SpectrographsProperties():
    """
    General properties of spectrographs.

    Attributes
    ----------
    data : pd.DataFrame
        Spectrograph properties read from file.
    dirout : str
        Output directory to save files and figures.
    
    Methods
    -------

    """
    def __init__(self, dirout='./', filprop=filprop):
        """
        Get spectrograph properties from file read with `read_spectrograph_properties` (default ./data/spectrograph/spectrograph_properties.csv)
        """
        self.data = read_spectrograph_properties(filprop=filprop)
        self.dirout = dirout

        if not os.path.exists(self.dirout): os.makedirs(self.dirout)


    def plot_wrange_line(self, ax=None, xunit='nm', xlabel=None, ylabel=None, title='', inst_in_plot=True, lisinst=None, yprop=None, cmap=None, va='bottom', lw=5):
        """
        Plot the wavelength range for each instrument.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            The axes on which to plot. If not provided, the current axes will be used.
        xunit : str, optional
            The unit of the x-axis. Default is 'nm'.
        xlabel : str, optional
            The label for the x-axis. If not provided, it will be generated based on the xunit.
        ylabel : str, optional
            The label for the y-axis.
        title : str, optional
            The title of the plot.
        inst_in_plot : bool, optional
            Whether to include the instrument acronym in the plot. Default is True. If False, it is recommended to add a legend.
        lisinst : list, optional
            A list of instrument names to plot. If not provided, all instruments will be plotted.
        yprop : str, optional
            The property to plot on the y-axis. If not provided, the y-axis will not be labeled.
        cmap : str, optional
            The colormap to use for coloring the lines. If not provided, lines will have a default color.
        va : str or list, optional
            The vertical alignment of the instrument acronym. If a list is provided, it should have the same length as lisinst. Default is 'bottom'.
        lw : float, optional
            The linewidth of the lines. Default is 5.

        Returns
        -------
        ax : matplotlib.axes.Axes
            The axes on which the plot is created.
        """
        if ax is None: ax = plt.gca()
        # Select wavelength unit
        wminstr, wmaxstr = f'wmin_{xunit}', f'wmax_{xunit}'
        # Select instruments
        if lisinst is None: df = self.data
        else: df = self.data.loc[lisinst]
        # Remove instrument with nan in properties to plot
        if yprop is None: masknan = [(np.isfinite(df.loc[i, wminstr])) and (np.isfinite(df.loc[i, wmaxstr])) for i in df.index]
        else: masknan = [(np.isfinite(df.loc[i, wminstr])) and (np.isfinite(df.loc[i, wmaxstr])) and (np.isfinite(df.loc[i, yprop])) for i in df.index]
        dfp = df[masknan]
        if isinstance(va, list): va = list(np.array(va)[masknan])

        # Colors
        if cmap is not None:
            cmap = mpl.colormaps.get_cmap(cmap)
            norm = mpl.colors.Normalize(vmin=0, vmax=len(dfp.index))

        # Plot
        for i, inst in enumerate(dfp.index):
            if yprop is None: yp = [-i, -i]
            else: yp = [dfp.loc[inst, yprop], dfp.loc[inst, yprop]]
            # Plot line
            if cmap is None:
                ax.plot([dfp.loc[inst, wminstr], dfp.loc[inst, wmaxstr]], yp, '-', label=dfp.loc[inst, 'inst_acronym_nice'], lw=lw, alpha=0.5)
            else:
                ax.plot([dfp.loc[inst, wminstr], dfp.loc[inst, wmaxstr]], yp, '-', label=dfp.loc[inst, 'inst_acronym_nice'], lw=lw, alpha=0.5, color=cmap(norm(i)))
            # Add instrument acronym text
            if isinstance(va, str): vap = va
            elif isinstance(va, list): vap = va[i]
            if inst_in_plot: ax.text(dfp.loc[inst, wminstr], yp[1], dfp.loc[inst, 'inst_acronym_nice'], va=vap, ha='left')

        # Make sure text within plot limits
        if yprop is None: 
            yp = [0.5, 0.5]  # largest yp 
            ax.plot([dfp.loc[inst, wminstr], dfp.loc[inst, wmaxstr]], yp, color='None')
        else:
            dy = (dfp[yprop].max() - dfp[yprop].min()) * 0.1
            ax.plot([dfp.loc[inst, wminstr], dfp.loc[inst, wmaxstr]], [dfp[yprop].max() + dy, dfp[yprop].max() + dy], color='None')

        # Style
        if xlabel is None: xlabel = wavelength_label(x=xunit)
        if (yprop is not None) and (ylabel is None): ylabel = dictprop_nice[yprop]
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
        if yprop is None:
            # Hide y ticks and labels
            ax.yaxis.set_tick_params(labelleft=False)
            ax.set_yticks([])
        elif (yprop == 'R') or (yprop == 'year_start') or (yprop == 'year_end') or (yprop == 'ndet') or (yprop == 'ndet_raw') or (yprop == 'nord'):
            # Force y-labels to be integers
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
        return ax


    def plot_wrange_rectangle(self, ax=None, xunit='nm', xlabel=None, ylabel=None, title='', inst_in_plot=True, lisinst=None, dyfrac=0.8, yprop=None, cmap=None, alpha=0.6, va='center'):
        """
        Plot the wavelength range for each instrument.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            The axes on which to plot. If not provided, the current axes will be used.
        xunit : str, optional
            The unit of the x-axis. Default is 'nm'.
        xlabel : str, optional
            The label for the x-axis. If not provided, it will be generated based on the xunit.
        ylabel : str, optional
            The label for the y-axis.
        title : str, optional
            The title of the plot.
        inst_in_plot : bool, optional
            Whether to include the instrument acronym in the plot. Default is True. If False, it is recommended to add a legend.
        lisinst : list, optional
            A list of instrument names to plot. If not provided, all instruments will be plotted.
        yprop : str, optional
            The property to plot on the y-axis. If not provided, the y-axis will not be labeled.
        cmap : str, optional
            The colormap to use for coloring the lines. If not provided, lines will have a default color.
        va : str or list, optional
            The vertical alignment of the instrument acronym. If a list is provided, it should have the same length as lisinst. Default is 'bottom'.
        lw : float, optional
            The linewidth of the lines. Default is 5.

        Returns
        -------
        ax : matplotlib.axes.Axes
            The axes on which the plot is created.
        """
        if ax is None: ax = plt.gca()
        # Select wavelength unit
        wminstr, wmaxstr = f'wmin_{xunit}', f'wmax_{xunit}'
        # Select instruments
        if lisinst is None: df = self.data
        else: df = self.data.loc[lisinst]
        # Remove instrument with nan in properties to plot
        if yprop is None: masknan = [(np.isfinite(df.loc[i, wminstr])) and (np.isfinite(df.loc[i, wmaxstr])) for i in df.index]
        else: masknan = [(np.isfinite(df.loc[i, wminstr])) and (np.isfinite(df.loc[i, wmaxstr])) and (np.isfinite(df.loc[i, yprop])) for i in df.index]
        dfp = df[masknan]
        if isinstance(va, list): va = list(np.array(va)[masknan])

        # Colors
        if cmap is not None:
            cmap = mpl.colormaps.get_cmap(cmap)
            norm = mpl.colors.Normalize(vmin=0, vmax=len(dfp.index))

        # Plot
        if yprop is None: dy = (len(dfp.index) - 0) / len(dfp.index) * dyfrac
        else: dy = np.abs(dfp[yprop].max() - dfp[yprop].min()) / len(dfp.index) * dyfrac
        for i, inst in enumerate(dfp.index):
            # Create a Rectangle patch
            x = dfp.loc[inst, wminstr]
            if yprop is None: y = -i
            else: y = dfp.loc[inst, yprop]
            w = dfp.loc[inst, wmaxstr] - dfp.loc[inst, wminstr]
            h = dy
            c = cmap(norm(i)) if cmap is not None else plt.rcParams['axes.prop_cycle'].by_key()['color'][i]
            rect = patches.Rectangle([x, y], w, h, linewidth=0, edgecolor='None', facecolor=c, alpha=alpha)
            # Add the patch to the Axes
            ax.add_patch(rect)
            # Add instrument acronym text
            if isinstance(va, str): vap = va
            elif isinstance(va, list): vap = va[i]
            if inst_in_plot: ax.text(dfp.loc[inst, wminstr], y+dy*0.5, dfp.loc[inst, 'inst_acronym_nice'], va=vap, ha='left')

        # Make sure text within plot limits
        if yprop is None:
            yp = [0.5, 0.5]  # largest yp 
            ax.plot([dfp.loc[inst, wminstr], dfp.loc[inst, wmaxstr]], yp, color='None')
        else:
            dy = (dfp[yprop].max() - dfp[yprop].min()) * 0.1
            ax.plot([dfp.loc[inst, wminstr], dfp.loc[inst, wmaxstr]], [dfp[yprop].max() + dy, dfp[yprop].max() + dy], color='None')

        # Style
        if xlabel is None: xlabel = wavelength_label(x=xunit)
        if (yprop is not None) and (ylabel is None): ylabel = dictprop_nice[yprop]
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
        if yprop is None:
            # Hide y ticks and labels
            ax.yaxis.set_tick_params(labelleft=False)
            ax.set_yticks([])
        elif (yprop == 'R') or (yprop == 'year_start') or (yprop == 'year_end') or (yprop == 'ndet') or (yprop == 'ndet_raw') or (yprop == 'nord'):
            # Force y-labels to be integers
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
        return ax


    def fig_wrange(self, filout=None, sh=False, sv=True, style='rectangle', **kwargs):
        """
        Plot the wavelength range of the spectrograph.

        Parameters
        ----------
        filout : str, optional
            The output file path to save the figure. If not provided, the figure will be saved in the directory specified by `self.dirout` with the filename `inst_wrange.pdf`.
        sh : bool, optional
            If True, display the figure using `plt.show()`. Default is False.
        sv : bool, optional
            If True, save the figure. Default is True.
        style : 'line', 'rectangle', optional
            Note: Thick lines (`lw`) can modify the wavelength range covered.
        **kwargs : dict, optional
            Additional keyword arguments to be passed to `self.plot_wrange()`.

        Returns
        -------
        None
        """
        fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)
        if style == 'line': ax = self.plot_wrange_line(ax=ax, **kwargs)
        elif style == 'rectangle': ax = self.plot_wrange_rectangle(ax=ax, **kwargs)
        if sh: plt.show()
        if sv:
            if filout is None: filout = f'{self.dirout}/inst_wrange.pdf'
            fig.savefig(filout)
        plt.close()
        return


class Spectrograph():
    """
    Spectrograph class

    Attributes
    ----------
    inst : str
        Instrument id. TODO: list of available instruments
    filprop : str
        File with spectrograph properties.
    ron:
        Readout noise per detector, red to blue [e-]
    """

    def __init__(self, inst, dirout='./', filprop=filprop, ordcut=True):
        """
        - Populate spectrograph properties from file read with `read_spectrograph_properties` (default ./data/spectrograph/spectrograph_properties.csv)
        - Get order data (wavelength range, order number, etc.) from file read with `read_ords` (default ./data/spectrograph/{inst}_orders.csv)

        Optional kwargs
        ordcut : bool
            Only applies to carmnir. If True, the orders are cut in two, for each detector, and treated as separate orders.
        """
        self.inst = inst
        self.filprop = filprop
        self.dirout = dirout

        if not os.path.exists(self.dirout): os.makedirs(self.dirout)

        # Get spectrograph properties from file
        dataspec = read_spectrograph_properties(filprop=self.filprop).loc[inst]
        self.R = dataspec['R']
        self.spectral_sampling_px = dataspec['spectral_sampling_px']
        self.pixel_ms = dataspec['pixel_ms']
        self.wmin_nm = dataspec['wmin_nm']
        self.wmax_nm = dataspec['wmax_nm']
        self.ndet = dataspec['ndet']
        self.ndet_raw = dataspec['ndet_raw']
        self.nord = dataspec['nord']
        self.ord_ref = dataspec['ord_ref']
        self.npix = dataspec['npix']
        self.nslice = dataspec['nslice']
        self.nslice_raw = dataspec['nslice_raw']
        self.ron = dataspec['ron']
        self.conversion_e_adu = dataspec['conversion_e_adu']
        self.position = dataspec['position']
        self.type = dataspec['type']
        self.inst_nice = dataspec['inst_nice']
        self.inst_acronym_nice = dataspec['inst_acronym_nice']
        self.tel = dataspec['tel']
        self.tel_nice = dataspec['tel_nice']
        self.tel_acronym_nice = dataspec['tel_acronym_nice']
        self.tel_diameter = dataspec['tel_diameter']
        self.observatory = dataspec['observatory']
        self.observatory_nice = dataspec['observatory_nice']
        self.observatory_acronym_nice = dataspec['observatory_acronym_nice']
        self.year_start = dataspec['year_start']
        self.year_end = dataspec['year_end']
        self.ref = dataspec['ref']
        self.notes = dataspec['notes']

        # Order data
        self.filords = get_filords(self.inst)
        self.dataord = read_ords(self.inst, self.filords, ordcut=ordcut)


    # def __repr__(self):
    #     return


    def __str__(self):
        # vars(self)
        return f'Instrument {self.inst} {self.inst_acronym_nice} ({self.inst_nice})'


    def print_properties_nice(self):
        for k in dictprop_nice.keys():
            print(f'{dictprop_nice[k]}: {vars(self)[k]}')
        return


    # -------------------------------------------------------------------------

    # Plotting

    def plot_ords_rectangle(self, ax=None, xunit='nm', nrows=2, olabel='ord_real', dyfrac=0.2, rowsep=0.2, xlabel=None, ylabel='', title='', legendlabel=None, legend=False, cmap=None, color='k', colorshading=True, va='bottom', ybase=0):
        """
        Plot orders

        Parameters
        ----------
        dyfrac : float
            Rectangle height fraction
        
        rowsep : float
            Extra separation between rows
        
        nrows : int
            1: All orders in the same row
            2: Alternate rows every each other order
            3: Alternate rows every three orders
            and so on
        olabel : 'ord_real', or 'ord_real_num'
            'ord_real': Real order number
            'ord_real_num': Real order number and order number from 0 to nord
        
        Returns
        -------
        ax : matplotlib.axes.Axes
            The axes on which the plot is created.
        """
        # ax = plt.gca() if ax is None else ax
        if ax is None: ax = plt.gca()
        if legendlabel is None: legendlabel = self.inst_acronym_nice

        dfp = self.dataord
        
        # Colors
        if cmap is not None:
            cmap = mpl.colormaps.get_cmap(cmap)
            norm = mpl.colors.Normalize(vmin=0, vmax=len(self.dataord.index))
        else:
            if colorshading:
                lisalpha = [1/nrows*nr for nr in np.arange(1, nrows+1, 1)]
                liscolor = [mpl.colors.to_rgba(color, a) for a in lisalpha]
            else:
                liscolor = [color for o in dfp.index]

        # Select wavelength unit
        wminstr, wmaxstr = f'wmin_{xunit}', f'wmax_{xunit}'

        # Plot
        for o in dfp.index:
            # Create a Rectangle patch
            x = dfp.loc[o, wminstr]
            y = ybase + o % nrows + (o % nrows) * rowsep
            w = dfp.loc[o, wmaxstr] - dfp.loc[o, wminstr]
            h = (len(dfp.index) - 0) / len(dfp.index) * dyfrac
            c = cmap(norm(o)) if cmap is not None else liscolor[o%nrows]
            rect = patches.Rectangle([x, y], w, h, linewidth=0, edgecolor='None', facecolor=c)
            # Add the patch to the Axes
            ax.add_patch(rect)
        
            # Add order label
            if isinstance(va, str): vap = va
            elif isinstance(va, list): vap = va[o]
            if olabel == 'ord_real':
                ostr = dfp.loc[o, 'ord_real']
                ax.text(x+w*0.5, y+h*0.9, ostr, va=vap, ha='center', fontsize='xx-small', color=c)
            elif olabel == 'ord_real_num':
                ostr = dfp.loc[o, 'ord_real']
                ax.text(x+w*0.5, y+h*0.9, ostr, va='bottom', ha='center', fontsize='xx-small', color=c)
                ostr = o
                ax.text(x+w*0.5, y-h*0.7, ostr, va='top', ha='center', fontsize='xx-small', color=c)

        # Fake plot to have correct limits
        ax.plot([dfp[wminstr].min(), dfp[wmaxstr].max()], [ybase, ybase + nrows + h*2 + rowsep*2], color='None')
        ax.plot([dfp[wminstr].min(), dfp[wmaxstr].max()], [ybase - h*2 - rowsep*2, ybase], color='None')
        
        ax.minorticks_on()
        # Hide y ticks and labels
        ax.yaxis.set_tick_params(labelleft=False)
        ax.set_yticks([])
        
        if xlabel is None: xlabel = wavelength_label(x=xunit)
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
        if legend: ax.legend()
        # xlim, ylim = ax.get_xlim(), ax.get_ylim()
        # dx, dy = xlim[1] - xlim[0], ylim[1] - ylim[0]
        # ax.set(xlim=(dfp[wminstr].min() - dx*0.02, dfp[wmaxstr].max() + dx*0.02), ylim=(ybase - dy*0.4, ybase+nrows-1 + dy*0.4))
        return ax


    def plot_ords_line(self, ax=None, xunit='nm', nrows=2, rowsep=0.2, olabel='ord_real', dyfrac=0.2, xlabel=None, ylabel='', title='', legendlabel=None, legend=False, cmap=None, color='k', va='bottom', ybase=0):
        """
        Plot orders

        Parameters
        ----------
        nrows : int
            1: All orders in the same row
            2: Alternate rows every each other order
            3: Alternate rows every three orders
            and so on
        olabel : 'ord_real', or 'ord_real_num'
            'ord_real': Real order number
            'ord_real_num': Real order number and order number from 0 to nord

        Returns
        -------
        ax : matplotlib.axes.Axes
            The axes on which the plot is created.
        """
        # ax = plt.gca() if ax is None else ax
        if ax is None: ax = plt.gca()
        if legendlabel is None: legendlabel = self.inst_acronym_nice

        dfp = self.dataord
        
        # Colors
        if cmap is not None:
            cmap = mpl.colormaps.get_cmap(cmap)
            norm = mpl.colors.Normalize(vmin=0, vmax=len(self.dataord.index))
        else:
            lisalpha = [1/nrows*nr for nr in np.arange(1, nrows+1, 1)]
            liscolor = [mpl.colors.to_rgba(color, a) for a in lisalpha]

        # Select wavelength unit
        wminstr, wmaxstr = f'wmin_{xunit}', f'wmax_{xunit}'

        # Plot
        for o in dfp.index:
            y = ybase + o % nrows + (o % nrows) * rowsep
            h = (len(dfp.index) - 0) / len(dfp.index) * dyfrac
            c = cmap(norm(o)) if cmap is not None else liscolor[o%nrows]

            # Wavelength range
            ax.annotate('', xy=(dfp.loc[o, wminstr], y + h), xytext=(dfp.loc[o, wmaxstr], y + h), arrowprops=dict(arrowstyle="|-|, widthA=0.25, widthB=0.25", lw=2, shrinkA=0, shrinkB=0, color=c), zorder=10)

            # Order label
            if olabel == 'ord_real':
                ax.text((dfp.loc[o, wminstr] + dfp.loc[o, wmaxstr]) / 2, y + h + 0.1, '{}'.format(dfp.loc[o, 'ord_real']), va='bottom', ha='center', fontsize='xx-small', color=c, bbox=dict(facecolor='white', edgecolor='None', alpha=0.3), zorder=1)
            elif olabel == 'ord_real_num':
                ax.text((dfp.loc[o, wminstr] + dfp.loc[o, wmaxstr]) / 2, y + h + 0.1, '{}'.format(dfp.loc[o, 'ord_real']), va='bottom', ha='center', fontsize='xx-small', color=c, bbox=dict(facecolor='white', edgecolor='None', alpha=0.3), zorder=1)
                ax.text((dfp.loc[o, wminstr] + dfp.loc[o, wmaxstr]) / 2, y + h - 0.1, '{}'.format(o), va='top', ha='center', fontsize='xx-small', color=c, bbox=dict(facecolor='white', edgecolor='None', alpha=0.3), zorder=1, rotation=0)

        # Fake plot to have correct limits
        ax.plot([dfp[wminstr].min(), dfp[wmaxstr].max()], [ybase, ybase + nrows + h*2 + rowsep*2], color='None')
        ax.plot([dfp[wminstr].min(), dfp[wmaxstr].max()], [ybase - h*2 - rowsep*2, ybase], color='None')
        
        ax.minorticks_on()
        # Hide y ticks and labels
        ax.yaxis.set_tick_params(labelleft=False)
        ax.set_yticks([])
        
        if xlabel is None: xlabel = wavelength_label(x=xunit)
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
        if legend: ax.legend()
        # xlim, ylim = ax.get_xlim(), ax.get_ylim()
        # dx, dy = xlim[1] - xlim[0], ylim[1] - ylim[0]
        # ax.set(xlim=(dfp[wminstr].min() - dx*0.02, dfp[wmaxstr].max() + dx*0.02), ylim=(ybase - dy*0.4, ybase+nrows-1 + dy*0.4))
        return ax


    def fig_ords(self, style='rectangle', filout=None, sh=False, sv=True, figsize=(16, 2), **kwargs):
        """
        """
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        if style == 'line': ax = self.plot_ords_line(ax=ax, **kwargs)
        elif style == 'rectangle': ax = self.plot_ords_rectangle(ax=ax, **kwargs)
        if sh: plt.show()
        if sv:
            if filout is None: filout = f'{self.dirout}/{self.inst}_ords.pdf'
            fig.savefig(filout)
        plt.close()
        return
