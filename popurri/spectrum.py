"""

"""
import os
from pathlib import Path
import sys

from astropy.io import fits
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

from .plotutils import wavelength_label
from . import spectrograph as spectrograph

dirhere = os.path.dirname(__file__)
dirdata = os.path.join(dirhere, '../data/')


###############################################################################


def read_spec_carmvis(filin):
    """
    Read CARACAL ouput.
    """
    dataspec = {}
    with fits.open(filin) as hdulist:
        dataspec['header'] = hdulist[0].header
        dataspec['w'] = hdulist['wave'].data
        dataspec['f'] = hdulist['spec'].data
        dataspec['fe'] = hdulist['sig'].data
        dataspec['c'] = hdulist['cont'].data
    # TODO check what is c exactly
    dataspec['nord'] = len(dataspec['w'])
    dataspec['ords'] = np.arange(0, dataspec['nord'], 1)
    dataspec['npix'] = len(dataspec['w'][0])
    return dataspec


def read_spec_carmnir(filin, ordcut=True, saveordnoncut=False):
    """
    NIR discontinuity at the center of the order: between pixel 2039 and pixel 2040. Divide orders in two.
    """
    dataspec = read_spec_carmvis(filin)

    if ordcut and saveordnoncut:
        for k in ['w', 'f', 'fe', 'c', 'ords']:
            dataspec[k+'_raw'] = dataspec[k].copy()
        for k in ['nord', 'npix']:
            dataspec[k+'_raw'] = dataspec[k]

    if ordcut:
        pixcut = int(dataspec['npix']/2)
        dataspec_new = {}
        for k in ['w', 'f', 'fe', 'c']:
            first = dataspec[k][:,:pixcut]
            second = dataspec[k][:,pixcut:]
            dataspec_new[k] = []
            for o in dataspec['ords']:
                dataspec_new[k].append(first[o])
                dataspec_new[k].append(second[o])
            dataspec[k] = np.array(dataspec_new[k])

        dataspec['nord'] = len(dataspec['w'])
        dataspec['ords'] = np.arange(0, dataspec['nord'], 1)
        dataspec['npix'] = len(dataspec['w'][0])
    
    # TODO This is the 2 detectors showing, think if can define thinks differently
    return dataspec


def read_spec_crires(filin):
    pass
    return


def read_spec_criresplus(filin):
    pass
    return


def read_spec_espresso(filin):
    pass
    return


# espresso_uhr11
# espresso_hr11
# espresso_hr21
# espresso_hr42
# espresso_mr42
# espresso_mr84
def read_spec_espresso_hr(filin):
    pass
    return


def read_spec_espresso_mr(filin):
    pass
    return


def read_spec_espresso_uhr(filin):
    pass
    return


def read_spec_expres(filin):
    pass
    return


def read_spec_harps(filin):
    pass
    return


def read_spec_harpsn(filin):
    pass
    return


def read_spec_maroonx(filin):
    pass
    return


def read_spec_neid(filin):
    pass
    return


def read_spectrum(filin, inst, **kwargs):
    """
    Read spectrum from file `filin`.

    Parameters
    ----------
    filin : str
        File with spectrum.
    inst : str
        Instrument keyword.
    
    ordcut=True
        carmnir
    saveordnoncut=False
        carmnir
    
    Returns
    -------
    dataspec : dict
        Dictionary with spectrum data.
    
    Notes
    -----

    dataspec general format (some keys might be empty/nan for some instruments):
    - 'w': Wavelength array
    - 'f': Flux array
    - 'fe': Flux error array
    - 'b': Blaze function array
    - 'm': Mask array
    - 'nord': Number of orders
    - 'ords': List of orders
    - 'header': Header dictionary
    - 'w1d': Wavelength array of the coadded 1D spectrum
    - 'f1d': Flux array of the coadded 1D spectrum
    - 'fe1d': Flux error array of the coadded 1D spectrum
    - 'b1d': Blaze function array of the coadded 1D spectrum
    - 'm1d': Mask array of the coadded 1D spectrum
    # TODO Change/generalise as more instruments are added
    """
    dictinst = {
        'carmvis': read_spec_carmvis,
        'carmnir': read_spec_carmnir,
        'crires': read_spec_crires,
        'criresplus': read_spec_criresplus,
        'espresso': read_spec_espresso,
        'espresso_hr': read_spec_espresso_hr,
        'express': read_spec_expres,
        'harps': read_spec_harps,
        'harpsn': read_spec_harpsn,
        'maroonx': read_spec_maroonx,
        'neid': read_spec_neid,
    }
    liskey = ['w', 'f', 'fe', 'b', 'm', 'nord', 'ords', 'npix', 'header', 'w1d', 'f1d', 'fe1d', 'b1d', 'm1d']

    # Check if file exists
    if not Path(filin).is_file(): sys.exit(f'File {filin} not found')

    # Remove extra kwargs
    if inst != 'carmnir':
        kwargs.pop('ordcut')
        kwargs.pop('saveordnoncut')

    # Read
    dataspec = dictinst[inst](filin, **kwargs)
    # Fill in missing keys
    for key in liskey:
        if key not in dataspec:
            dataspec[key] = None

    return dataspec


def read_spectra(lisfilin, inst, returnclass=False, dirout='./', **kwargs):
    """
    returnclass : bool
        Return a list of Spectrum objects if True, or a list of dictionaries if False.
    kwargs
    dirout : Needed if returnclass=True
    ordcut : bool
        carmnir
    saveordnoncut : bool
        carmnir
    """
    lisdataspec = []
    m = np.ones(len(lisfilin), dtype=bool)
    for i, filin in enumerate(lisfilin):
        try:
            if returnclass:
                lisdataspec.append(Spectrum(filin, inst, dirout, **kwargs))
            else:
                lisdataspec.append(read_spectrum(filin, inst, **kwargs))
        except:
            print(f'Error reading {filin}')
            m[i] = False
        # # TESTING
        # if returnclass:
        #     print(f'Reading {filin}')
        #     lisdataspec.append(Spectrum(filin, inst, dirout, **kwargs))
        # else:
        #     lisdataspec.append(read_spectrum(filin, inst, **kwargs))

    lisfilin = np.array(lisfilin)[m]
    return lisdataspec, lisfilin


def read_header(filin, ext=0):
    """
    Read FITS header. Can choose which extension with `ext` (default 0).
    """
    with fits.open(filin) as hdulist:
        header = hdulist[ext].header
    return header


def read_header_kw_table():
    """
    Read table with parameters keyword in this code and FITS header keywords from file 'header_kws.csv' in directory `dirdata`.
    """
    filtable = os.path.join(dirdata, 'spectrograph/header_kws.csv')
    df = pd.read_csv(filtable, comment='#').set_index('param')
    return df


def get_params_header(filin, inst, headertable=None, ext=0):
    """
    Get general common parameters from header.

    Table with parameters keyword in this code and FITS header keywords in file 'header_kws.csv' in directory `dirdata`.

    filin : str or astropy header object (astropy.io.fits.header.Header)
        FITS file from which to read the header or astropy header object (faster if header already loaded before). If it is a header object, the parameter `ext` is not used.
    """
    # Read header
    if isinstance(filin, str): header = read_header(filin, ext=ext)
    elif isinstance(filin, fits.header.Header): header = filin

    # Get reference table
    if headertable is None: headertable = read_header_kw_table()
    headertable = headertable[inst]

    # Get keyword values
    data = {}
    for p in headertable.keys():
        if headertable[p] is not np.nan:
            data[p] = header[headertable[p]]
        else:
            data[p] = np.nan
    return data
            

def get_header_snr(filin, inst, ext=0):
    """
    """
    # Read header
    if isinstance(filin, str): header = read_header(filin, ext=ext)
    elif isinstance(filin, fits.header.Header): header = filin

    # S/N global keyword
    dictpattern = {
        'carmvis': '*CARACAL FOX SNR*',
        'carmnir': '*CARACAL FOX SNR*',
        # TODO add all instruments
        # TODO carmnir SNR 55 (ordcut) orders, not clear which is which
    }
    # kws = ['HIERARCH CARACAL FOX SNR {:d}'.format(o) for o in ords]

    pattern = dictpattern[inst]
    lish = header[pattern]  # astropy.io.fits.header.Header object
    vals = list(lish.values())
    kws = list(lish.keys())
    # Change kws to order number
    kws = np.arange(0, len(vals), 1)  # NOTE this assumes vals are in order from blue to red order
    dictsnr = dict(zip(kws, vals))  # dict object
    return dictsnr


###############################################################################

class Spectrum():
    """
    Spectrum class

    dataspec?

    dataheader
    dataord
    """

    def __init__(self, filin, inst, dirout='./', obj=None, tag=None, ordcut=True, saveordnoncut=False, headertable=None):
        """

        ordcut : bool
            carmnir
        saveordnoncut : bool
            carmnir
        """
        self.filin = filin
        self.filname = os.path.basename(os.path.splitext(filin)[0])
        self.inst = inst
        self.obj = obj
        self.tag = tag
        self.dirout = dirout
        if not os.path.exists(self.dirout): os.makedirs(self.dirout)

        # Read spectrum
        dataspec = read_spectrum(self.filin, self.inst, ordcut=ordcut, saveordnoncut=saveordnoncut)
        # Get main (array) data as attributes
        # TODO Change and keep as dataspec dictionary?
        self.w = dataspec['w']
        self.f = dataspec['f']
        self.fe = dataspec['fe']
        self.b = dataspec['b']
        self.m = dataspec['m']
        self.nord = dataspec['nord']
        self.ords = dataspec['ords']
        self.npix = dataspec['npix']
        self.header = dataspec['header']
        self.w1d = dataspec['w1d']
        self.f1d = dataspec['f1d']
        self.fe1d = dataspec['fe1d']
        self.b1d = dataspec['b1d']
        self.m1d = dataspec['m1d']

        # Add Spectrograph object
        self.Spectrograph = spectrograph.Spectrograph(self.inst, dirout=self.dirout, ordcut=ordcut)
        # Add "real" orders from Spectrograph object
        self.ords_real = self.Spectrograph.dataords['ord_real'].values

        # carmnir: cut orders in two. Add extra attributes if needed
        if (self.inst == 'carmnir') and (ordcut is True) and (saveordnoncut is True):
            self.w_raw = dataspec['w_raw']
            self.f_raw = dataspec['f_raw']
            self.fe_raw = dataspec['fe_raw']
            self.c_raw = dataspec['c_raw']
            self.nord_raw = dataspec['nord_raw']
            self.ords_raw = dataspec['ords_raw']
            self.npix_raw = dataspec['npix_raw']

        # Add object from header
        if self.obj is None:
            self.obj = self.header['OBJECT']  # TODO check if same for all instruments

        # Get parameters from header in `dataheader` attribute
        self.dataheader = get_params_header(self.header, self.inst, headertable=headertable)

        # Get S/N per order
        # TODO
        self.dataord = pd.DataFrame({
            'snr': get_header_snr(self.filin, self.inst, ext=0)
        })



    # def __repr__(self):
    #     return


    def __str__(self):
        # vars(self)
        return f''
    

    def plot_spectrum(self, ax=None, ords=None, wmin=None, wmax=None, legend=False, legendloc=None, xunit='A', xlabel=None, ylabel='Flux', title='', lw=1, linestyle='-', alpha=1, zorder=0):
        """Plot spectrum flux vs wavelength (or pixel), for the orders in `ords`.
        
        Parameters
        ----------
        ords: list-like, optional
            Orders to plot. If None, plot all orders.
        wmin, wmax: floats, optional
            Wavelength range to plot. Must be in the same units as `w` (and `xunit`).
        TODO Add support for color from colormap
        """
        if ax is None: ax = plt.gca()
        if ords is None: ords = self.ords
        if np.issubdtype(type(ords), np.integer): ords = [ords]  # make sure it is a list

        # Wavelength range
        if wmin is None: wmin = np.nanmin(self.w[ords])
        if wmax is None: wmax = np.nanmax(self.w[ords])
        mp = (self.w >= wmin) & (self.w <= wmax)
        
        # Plot
        for o in ords:
            mpo = mp[o]
            if not any(mpo): continue  # skip if all pixels masked
            ax.plot(self.w[o][mpo], self.f[o][mpo], lw=lw, linestyle=linestyle, alpha=alpha, zorder=zorder, label=f'{o}')
        if legend: ax.legend(loc=legendloc)
        if xlabel is None: xlabel = wavelength_label(x=xunit)
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title, xlim=(wmin, wmax))
        ax.minorticks_on()
        return ax
    

    def fig_spectrum(self, filout=None, sh=False, sv=True, figsize=(16, 4), **kwargs):
        """
        """
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        ax = self.plot_spectrum(ax=ax, **kwargs)
        if sh: plt.show()
        if sv:
            if filout is None: filout = f'{self.dirout}/{self.filname}_spec.pdf'
            fig.savefig(filout)
        plt.close()
        return
    

    def plot_spectrum_map(self, ax=None, ords=None, pixs=None, xlabel='Pixel', ylabel='Order', cblabel='Flux', title='', ordtype='zero', vmin=None, vmax=None, cmap='viridis'):
        """Plot spectrum flux matrix, all orders (orders in y-axis, pixel in x-axis).
        ords: list-like, optional
            Orders to plot (in zero-based notation). If None, plot all orders. Cannot be a single order.
        pixs: list-like, optional
            Pixels to plot. If None, plot all pixels. Cannot be a single pixel.
        ordtype : 'real' or 'zero', default 'real'
            If 'real', the real order number is plotted. If 'zero', the zero-based order number is plotted.
        """
        # TODO imshow doesn't really get the centre of the axis correcty. Try pcolormesh or matshow instead
        if ax is None: ax = plt.gca()
        # Orders to plot
        if ords is None: mords = np.ones_like(self.ords, dtype=bool)
        else:
            mords = np.zeros_like(self.ords, dtype=bool)
            mords[ords] = True
        # Order labels
        if ordtype == 'real': ords = self.ords_real[mords]
        elif ordtype == 'zero': ords = self.ords[mords]
        # Pixels to plot
        if pixs is None:
            pixs = np.arange(0, self.npix, 1)
            mpixs = np.ones_like(self.f[0], dtype=bool)
        else:
            mpixs = np.zeros_like(self.f[0], dtype=bool)
            mpixs[pixs] = True
        extent = [pixs[0], pixs[-1], ords[0], ords[-1]]
        # Flux limits
        if (vmin is not None) and (vmax is not None): extend = 'both'
        elif (vmin is not None) and (vmax is None): extend = 'min'
        elif (vmin is None) and (vmax is not None): extend = 'max'
        else: extend = None
        im = ax.imshow(self.f[mords,:][:,mpixs], origin='lower', interpolation='none', aspect='auto', extent=extent, vmin=vmin, vmax=vmax, cmap=cmap, rasterized=True)
        cbar = plt.colorbar(im, extend=extend, label=cblabel)
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
        # Force y-labels to be integers
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))        
        return ax


    def fig_spectrum_map(self, filout=None, sh=False, sv=True, figsize=(8, 6), **kwargs):
        """
        """
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        ax = self.plot_spectrum_map(ax=ax, **kwargs)
        if sh: plt.show()
        if sv:
            if filout is None: filout = f'{self.dirout}/{self.filname}_specmap.pdf'
            fig.savefig(filout)
        plt.close()
        return


    def plot_dataords(y, z=None):
        pass
        # TODO!


class Spectra():
    """
    Collection of (time series) spectra of the same instrument.
    """

    def __init__(self, lisfilin, inst, dirout='./', obj=None, tag=None, ordcut=True, saveordnoncut=False, deleteindividual=False, ext=0):
        """
        tag : str, optional
            Tag to add to the object name (e.g. type of reduction, sky subtracted...).

        ext : int, optional
            Header extension to read.

        Accessing the properties of Spectra
        ```
        # `lisfilin` is a list of 4 spectra with 61 orders and 4096 pixels each
        
        # Acces the flux `f` of all spectra
        >>> spectra = Spectra(lisfilin, inst)
        >>> spectra.f.shape
        (4, 61, 4096)
        
        # Access the flux `f` of order 4 for all spectra:
        >>> spectra.f[:,4,:].shape
        (4, 4096)

        # Access the flux `f` of orders 4 to 8 and pixels 1500 to 1600 for the first 3 observations:
        >>> spectra.f[:3,4:9,1500:1601].shape
        (3, 5, 101)
        ```
        """
        self.lisfilin = lisfilin
        self.lisfilin_initial = lisfilin
        self.inst = inst
        self.obj = obj
        self.tag = tag
        self.dirout = dirout

        if not os.path.exists(self.dirout): os.makedirs(self.dirout)

        headertable = read_header_kw_table()

        # Read files into Spectrum objects
        self.lisspec, self.lisfilin = read_spectra(self.lisfilin_initial, self.inst, returnclass=True, dirout=dirout, ordcut=ordcut, saveordnoncut=saveordnoncut, headertable=headertable)

        # Rearrange individial data into property for all spectra
        self.filname = np.array([spec.filname for spec in self.lisspec])
        self.w = np.array([spec.w for spec in self.lisspec])
        self.f = np.array([spec.f for spec in self.lisspec])
        self.fe = np.array([spec.fe for spec in self.lisspec])
        self.b = np.array([spec.b for spec in self.lisspec])
        self.m = np.array([spec.m for spec in self.lisspec])
        self.header = [spec.header for spec in self.lisspec]
        self.w1d = np.array([spec.w1d for spec in self.lisspec])
        self.f1d = np.array([spec.f1d for spec in self.lisspec])
        self.fe1d = np.array([spec.fe1d for spec in self.lisspec])
        self.b1d = np.array([spec.b1d for spec in self.lisspec])
        self.m1d = np.array([spec.m1d for spec in self.lisspec])

        # For general properties shared by all spectra, add only the value from the first spectrum
        self.nord = self.lisspec[0].nord
        self.ords = self.lisspec[0].ords
        self.ords_real = self.lisspec[0].ords_real
        self.npix = self.lisspec[0].npix
        if self.obj is None: self.obj = self.lisspec[0].obj
        if self.tag is None: self.tag = self.lisspec[0].tag
        
        # carmnir: cut orders in two. Add extra attributes if needed
        if (self.inst == 'carmnir') and (ordcut is True) and (saveordnoncut is True):
            self.w_raw = np.array([spec.w_raw for spec in self.w_raw])
            self.f_raw = np.array([spec.f_raw for spec in self.f_raw])
            self.fe_raw = np.array([spec.fe_raw for spec in self.fe_raw])
            self.c_raw = np.array([spec.c_raw for spec in self.c_raw])
            self.nord_raw = np.array([spec.nord_raw for spec in self.nord_raw])
            self.ords_raw = np.array([spec.ords_raw for spec in self.ords_raw])
            self.npix_raw = np.array([spec.npix_raw for spec in self.npix_raw])

        # Get parameters from header in `dataheader` attribute
        # Transform list of dictionaries into pandas dataframe
        lisdataheader = [sp.dataheader for sp in self.lisspec]
        self.dataheader = pd.DataFrame(lisdataheader, index=self.filname)
        del lisdataheader

        # Delete individual Spectrum objects to save memory
        if deleteindividual: del self.lisspec

    # def __repr__(self):
    #     return


    def __str__(self):
        # vars(self)
        return f''
    

    def plot_spectra(self, ax=None, ords=None, lisspec=None, wmin=None, wmax=None, legendlabel=None, legendloc=None, xunit='A', xlabel=None, ylabel='Flux', title='', cprop=None, cprop_all=False, cbar=None, cbarlabel=None, cmap=None, lw=1, linestyle0='-', linestyle1='-', alpha0=1, alpha1=0.7):
        """
        Plot spectra flux vs wavelength, orders in `ords`.

        Plot different observations with different colors, and optionally alternating alpha for different orders.

        Parameters
        ----------
        ords : list-like, optional
            List of orders to plot. If `None` (default), all orders are plotted.
        
        lisspec : list-like, optional
            List of observations (the number) to plot. Default is None, which plots all spectra in `lisspec`.

        cprop : str, optional
            Color based on property of the spectrum (e.g. S/N, airmass, BJD...). Uses `cmap` to define the color map. If `cmap` is None, use colormap viridis.

        cprop_all : bool, optional
            If True, get the colorbar limits from all spectra (even the ones not plotted). If False, adapt the colorbar limits to the plotted spectra.
            
        cbar :
            Add color bar. Defaults to True if cprop is not None, else defaults to False.

        wmin, wmax: floats, optional
            Wavelength range to plot. Must be in the same units as `w` (and `xunit`).

        linestyel0, linestyle1 : str, optional
            For alternating orders
        alpha0, alpha1 : float, optional
            For alternating orders
            
        legendlabel :
            None
            'obs'
        """
        if ax is None: ax = plt.gca()
        if lisspec is None: lisspec = np.arange(len(self.lisspec))
        if np.issubdtype(type(lisspec), np.integer): lisspec = [lisspec]  # make sure it is a list
        if ords is None: ords = self.ords
        if np.issubdtype(type(ords), np.integer): ords = [ords]  # make sure it is a list
        
        # Wavelength range
        if wmin is None: wmin = np.nanmin(self.w[lisspec][:,ords])
        if wmax is None: wmax = np.nanmax(self.w[lisspec][:,ords])
        mp = (self.w >= wmin) & (self.w <= wmax)

        # if cprop is not None:
        #     if cmap is None: cmap = 'viridis'
        #     cmap = mpl.colormaps.get_cmap(cmap)
        #     norm = mpl.colors.Normalize(vmin=np.nanmin(cprop), vmax=np.nanmax(cprop))

        # Colors
        if cmap is not None:
            cmap = mpl.colormaps.get_cmap(cmap)
            if cprop is not None:  # use property to define color
                if cprop_all:
                    norm = mpl.colors.Normalize(vmin=np.nanmin(cprop), vmax=np.nanmax(cprop))
                else:
                    norm = mpl.colors.Normalize(vmin=np.nanmin(cprop[lisspec]), vmax=np.nanmax(cprop[lisspec]))
            else:  # use spectrum number to define color
                norm = mpl.colors.Normalize(vmin=0, vmax=len(lisspec))

        # Color bar
        if (cbar is None) and (cprop is not None): cbar = True
        if cbar is None: cbar = False

        # Plot
        for i in lisspec:
            spec = self.lisspec[i]
            if cprop is not None:
                c = cmap(norm(cprop[i]))
            elif cmap is not None:
                c = cmap(norm(i))
            else:
                c = plt.rcParams['axes.prop_cycle'].by_key()['color'][i]
            # c = cmap(norm(i)) if cmap is not None else plt.rcParams['axes.prop_cycle'].by_key()['color'][i]
            # Loop over orders
            for o in ords:
                mpo = mp[i][o]
                if not any(mpo): continue  # skip if all pixels masked
                a = alpha0 if o % 2 == 1 else alpha1
                linestyle = linestyle0 if o % 2 == 1 else linestyle1
                label = f'{spec.filname}' if o == ords[0] else None
                ax.plot(spec.w[o][mpo], spec.f[o][mpo], c=c, alpha=a, lw=lw, linestyle=linestyle, label=label,)
        if cbar:
            cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label=cbarlabel)
            cb.minorticks_on()
        if legendlabel is not None: ax.legend(loc=legendloc)
        if xlabel is None: xlabel = wavelength_label(x=xunit)
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title, xlim=(wmin, wmax))
        ax.minorticks_on()
        return ax


    def fig_spectra(self, filout=None, sh=False, sv=True, figsize=(16, 4), **kwargs):
        """
        """
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        ax = self.plot_spectra(ax=ax, **kwargs)
        if sh: plt.show()
        if sv:
            if filout is None:
                if self.tag is not None: 
                    filout = f'{self.dirout}/{self.tag}_spectra.pdf'
                else:
                    filout = f'{self.dirout}/{self.obj}_spec.pdf'
            fig.savefig(filout)
        plt.close()
        return
    

    def plot_spectra_map(self, o=0):
        """
        Plot spectra flux matrix, order in `o`.
        """
        pass
        return






