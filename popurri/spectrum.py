"""
`spectrum` module

Methods to read and plot spectra from different instruments.
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
dirdata = os.path.join(dirhere, './data/')


###############################################################################

# Functions to read spectra and extra data
# ----------------------------------------

# CARMENES

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
    dataspec['pix'] = np.array([np.arange(0, dataspec['npix'], 1)] * len(dataspec['ords']))
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


# =============================================================================

# CRIRES

def read_spec_crires(filin):
    pass
    return

# =============================================================================

def read_spec_criresplus(filin):
    pass
    return


# =============================================================================

# ESPRESSO

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


# =============================================================================

# EXPRES

def read_spec_expres(filin):
    pass
    return


# =============================================================================

# HARPS

def wpolycoeff(filin, kwinst='HIERARCH ESO', nord=72):
    """From e2ds header FITS file get the polynomial coefficients necessary to transform from pixels to wavelength.

    Parameters
    ----------
    filin : str or astropy header object (astropy.io.fits.header.Header)
        FITS file from which to read the header or astropy header object (faster if header already loaded before).
    """

    # Read header
    if isinstance(filin, str):
        with fits.open(filin) as hdulist:
            dataspec['header'] = hdulist[0].header
    elif isinstance(filin, fits.header.Header):
        header = filin

    ords = np.arange(0, nord, 1)

    # Get polynomial coefficients (4 per order) from header
    #  Polynomial coefficient numbers, 4 per order
    polydeg = header[f'{kwinst} DRS CAL TH DEG LL']
    coeffnum = [np.arange(0+o*(polydeg+1), (polydeg+1)+o*(polydeg+1), 1) for o in ords]
    #  Polynomial coefficient values, 4 per order
    coeff = [[header[kwinst + 'DRS CAL TH COEFF LL{:d}'.format(j)] for j in coeffnum[o]] for o in ords]

    return coeff


def pix2wave(x, coeff):
    """Convert pixel to wavelength using the coefficients from the e2ds header, for a single echelle order.

    Parameters
    ----------
    x : 1d array
        Pixels.
    """
    w = coeff[0] + coeff[1] * x + coeff[2] * x**2 + coeff[3] * x**3
    return w


def pix2wave_echelle(x, coeff):
    nord = len(coeff)
    ords = np.arange(0, nord, 1)
    w = np.array([pix2wave(x, coeff[o]) for o in ords])
    return w



def read_spec_harps(filin, kwinst='HIERARCH ESO', nord=72, readblaze=True, dirblaze=None, filblaze=None):
    """
    Read e2ds reduced spectrum flux and wavelength, and optionally the blaze.

    The wavelength data is obtained from the header keywords in `filin` using `pix2wave_echelle`, `pix2wave`, and `wpolycoeff`.

    The blaze is only read if `readblaze` is True (default).
    The directory containing the blaze files by default is the same as the directory where the e2ds file is (`filin`), but can be changed with `dirblaze`.
    The blaze file by default is obtained from the header of the e2ds file `filin`: 'HIERARCH TNG DRS BLAZE FILE', but can be changed with `filblaze`.

    Parameters
    ----------
    filin : str
        Reduced e2ds file.
    readblaze : bool, default True
        Whether to read the blaze or not. If False, the returned blaze `b` is an array full of ones of the same shape as the spectrum.
    dirblaze : str, default None
        Directory where the blaze file is. If None (default), it is assumed that it is in the same directory as the spectrum `filin`.
    filblaze : str, default None
        Blaze file. Use if want to obtain the blaze from a file different than the one specified in the header keyword 'HIERARCH ESO DRS BLAZE FILE'.
    """
    dataspec = {}

    # Flux
    with fits.open(filin) as hdulist:
        dataspec['header'] = hdulist[0].header
        dataspec['f'] = hdulist[0].data

    # Wavelength
    # nord = len(dataspec['f'])
    npix = len(dataspec['f'][0])
    pix = np.arange(0, npix, 1)
    coeff = wpolycoeff(dataspec['header'], kwinst=kwinst, nord=nord)
    dataspec['w'] = pix2wave_echelle(x, coeff)

    # Blaze function
    if readblaze:
        # If no blaze directory specified, assume is the same as the data file
        if dirblaze is None: dirblaze = os.path.dirname(filin)
        # If no blaze file specified, get the corresponding one from the ehader
        if filblaze is None: filblaze = dataspec['header'][f'{kwinst} DRS BLAZE FILE']
        with fits.open(os.path.join(dirblaze, filblaze)) as hdulist:
            dataspec['b'] = hdulist[0].data
    else:
        dataspec['b'] = np.ones_like(dataspec['f'])        

    dataspec['nord'] = nord
    dataspec['ords'] = np.arange(0, dataspec['nord'], 1)
    dataspec['npix'] = npix
    dataspec['pix'] = np.array([np.arange(0, dataspec['npix'], 1)] * len(dataspec['ords']))
    return dataspec


# HARPS-N

def read_spec_harpsn(filin, kwinst='HIERARCH TNG', nord=69, readblaze=True, dirblaze=None, filblaze=None):
    """
    Read HARPS-N e2ds reduced spectrum flux and wavelength, and optionally the blaze.
    """
    dataspec = read_spec_harps(filin, kwinst=kwinst, nord=nord, readblaze=readblaze, dirblaze=dirblaze, filblaze=filblaze)
    return dataspec


# =============================================================================

# MAROON-X

def read_spec_maroonx(filin):
    pass
    return



# =============================================================================

# NEID

def read_spec_neid(filin):
    pass
    return


# dict_read_spectrum = {}
# def register_inst_read_spectrum(func):
#     """Decorator to register a function to read a spectrum for a specific instrument with the general function `read_spectrum`."""
#     dict_read_spectrum[func.__name__] = func
#     return func


# =============================================================================

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
    dict_read_spectrum = {
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
    liskey = ['w', 'f', 'fe', 'b', 'c', 'm', 'pix', 'nord', 'ords', 'npix', 'header', 'w1d', 'f1d', 'fe1d', 'b1d', 'm1d']

    # Check if file exists
    if not Path(filin).is_file(): sys.exit(f'File {filin} not found')

    # Remove extra kwargs
    if inst != 'carmnir':
        kwargs.pop('ordcut')
        kwargs.pop('saveordnoncut')
    if (inst != 'harps') and (inst != 'harpsn'):
        kwargs.pop('readblaze')
        kwargs.pop('dirblaze')
        kwargs.pop('filblaze')

    # Read
    dataspec = dict_read_spectrum[inst](filin, **kwargs)

    # Fill in missing keys
    for key in liskey:
        if key not in dataspec:
            dataspec[key] = None

    return dataspec


def read_spectra(lisfil, inst, returnclass=False, dirout='./', **kwargs):
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
    # Check that lisfil is a list or array, and not a string with the name of a single file
    if isinstance(lisfil, str):
        lisfil = [lisfil]

    lisdataspec = []
    m = np.ones(len(lisfil), dtype=bool)
    for i, filin in enumerate(lisfil):
        try:
            if returnclass:
                lisdataspec.append(Spectrum(filin, inst, dirout, **kwargs))
            else:
                lisdataspec.append(read_spectrum(filin, inst, **kwargs))
        except FileNotFoundError:
            print(f'File {filin} not found, removing from list')
            m[i] = False
    lisfil = np.array(lisfil)[m]
    return lisdataspec, lisfil


def read_header(filin, ext=0):
    """
    Read FITS header. Can choose which extension with `ext` (default 0).
    """
    with fits.open(filin) as hdulist:
        header = hdulist[ext].header
    return header


filheader_kw_table = os.path.join(dirdata, 'spectrograph/header_kws.csv')

def read_header_kw_table(filtable=filheader_kw_table, verbose=False):
    """
    Read table that contains general parameters keywords used in this code and their corresponding FITS header keyword for different instruments. The table is in file 'header_kws.csv' in directory `dirdata`.

    Short example of the table structure:
    ```
    |Param   |carmvis              |harps                |
    |--------|---------------------|---------------------|
    |airmass |AIRMASS              |AIRMASS              |
    |berv    |HIERARCH CARACAL BERV|HIERARCH ESO DRS BERV|
    ```
    """
    df = pd.read_csv(filtable, comment='#').set_index('param')
    if verbose: print('General parameters from FITS header:', df.index)
    return df


def get_params_header(filin, inst, headertable=None, ext=0):
    """
    Get general common parameters from header.

    The general parameters are taken from the table 'header_kws.csv' in directory `dirdata`, read with `read_header_kw_table`.

    filin : str or astropy header object (astropy.io.fits.header.Header)
        FITS file from which to read the header or astropy header object (faster if header already loaded before). If it is a header object, the parameter `ext` is not used.
    headertable : pandas DataFrame, optional
        Table previously read with `read_header_kw_table`. If None, it is read now.
    """
    # Get reference table
    if headertable is None: headertable = read_header_kw_table()
    headertable = headertable[inst]

    # Read header
    if isinstance(filin, str): header = read_header(filin, ext=ext)
    elif isinstance(filin, fits.header.Header): header = filin

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
    Keywords are the order number (starting at 0 for the bluest order).
    """
    # Read header
    if isinstance(filin, str): header = read_header(filin, ext=ext)
    elif isinstance(filin, fits.header.Header): header = filin

    # S/N global keyword
    dictpattern = {
        'carmvis': '*CARACAL FOX SNR*',
        'carmnir': '*CARACAL FOX SNR*',
        'harps': '*ESO DRS SPE EXT SN*',
        'harpsN': '*TNG DRS SPE EXT SN*',
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


# Other utils
# -----------

# Wavelength vacuum/air

def wvac2air(w):
    """Transform vacuum wavelength to air wavelength.
    Formula from: Ciddor 1996, Applied Optics 62, 958.

    Parameters
    ----------
    w : float or array-like
        Vacuum wavelength to be transformed to air, in A. If array-like, w sorted in increasing or decreasing order.
    """
    scalar = False
    if isinstance(w, (int, float)):
        w = [w]
        scalar = True
    w = np.array([w])
    wair = w.copy()

    mask = w > 2000.  # Modify only wavelength above 2000 A

    s2 = (1e4/w[mask])**2
    f = 1.+0.05792105/(238.0185-s2)+0.00167917/(57.362-s2)
    wair[mask] = w[mask]/f
    return wair[0][0] if scalar else wair[0]


def wair2vac(w):
    """Transform air wavelength to vacuum wavelength.
    Formula from: Ciddor 1996, Applied Optics 62, 958.

    w : float or array-like
        Air wavelength to be transformed to vacuum, in A. If array-like, w sorted in increasing or decreasing order.
    """
    scalar = False
    if isinstance(w, (int, float)):
        w = [w]
        scalar = True
    w = np.array([w])
    wvac = w.copy()

    mask = w > 2000.  # Modify only wavelength above 2000 A

    s2 = (1e4/w[mask])**2
    f = 1.+0.05792105/(238.0185-s2)+0.00167917/(57.362-s2)
    wvac[mask] = w[mask]*f
    return wvac[0][0] if scalar else wvac[0]


# Doppler shift

def dopplershift(x, v, rel=True):
    """
    x : float
        Wavelength
    v : float
        Velocity, m/s
    """
    C_MS = 2.99792458*1.e8  # Light speed [m/s]
    if rel: a = np.sqrt((1 + v / C_MS) / (1 - v / C_MS))
    else: a = (1 - v / C_MS)
    xprime = x * a
    return xprime



###############################################################################


# Spectrum observations classes
# -----------------------------

class Spectrum():
    """
    Spectrum class

    Data attributes of the `Spectrum` class:
    - `dataspec`: Dictionary with spectrum data (wavelength, flux, blaze...).
    - `dataheader`: Dictionary with general parameters from the header.
    - `dataord`: Pandas DataFrame with per order data such as the S/N.
    """

    def __init__(self, filin, inst, dirout='./', obj=None, tag='', headertable=None, ord_ref=None,
    # carmnir parameters
    ordcut=True, saveordnoncut=False,
    # harps and harpsn parameters
    readblaze=True, dirblaze=None, filblaze=None,
    ):
        """

        dataspec
        header
        dataheader
        dataords
        
        carmnir parameters
        ------------------
        ordcut : bool
        saveordnoncut : bool
        
        harps and harpsn parameters
        ---------------------------
        readblaze
        dirblaze
        filblaze
        """
        self.filin = filin
        self.filname = os.path.basename(os.path.splitext(filin)[0])
        self.inst = inst
        self.obj = obj
        self.tag = tag
        self.dirout = dirout
        if not os.path.exists(self.dirout): os.makedirs(self.dirout)

        # Read spectrum in `dataspec` attribute (dictionary)
        self.dataspec = read_spectrum(self.filin, self.inst, ordcut=ordcut, saveordnoncut=saveordnoncut, readblaze=readblaze, dirblaze=dirblaze, filblaze=filblaze)

        # Add Spectrograph object
        self.Spectrograph = spectrograph.Spectrograph(self.inst, dirout=self.dirout, ordcut=ordcut)
        # Add "real" orders from Spectrograph object
        self.ords_real = self.Spectrograph.dataord['ord_real'].values
        # Add order reference if ord_ref is None
        self.ord_ref = self.Spectrograph.ord_ref if ord_ref is None else ord_ref
        # Add pixel size in velocity [m/s]
        self.pixel_ms = self.Spectrograph.pixel_ms

        # Add general parameters from `dataspec` attribute (and remove from self.dataspec)
        self.nord = self.dataspec.pop('nord', None)
        self.ords = self.dataspec.pop('ords', None)
        self.npix = self.dataspec.pop('npix', None)
        self.header = self.dataspec.pop('header', None)

        # Add object from header
        if self.obj is None:
            self.obj = self.header['OBJECT']  # TODO check if same for all instruments

        # Get parameters from header in `dataheader` attribute (dictionary)
        self.dataheader = get_params_header(self.header, self.inst, headertable=headertable)

        # Get per order parameters in `dataord` attribute (pandas DataFrame)
        # TODO
        self.dataord = pd.DataFrame({
            'snr': get_header_snr(self.filin, self.inst, ext=0)
        })

        # Add per order parameters from `dataord` to `dataheader`. Keywords will be '{parameter}o{ord}' where parameter is the column in `dataord` (e.g. 'snr') and ord is the order number.
        for col in self.dataord.columns:
            for o in self.dataord.index:
                self.dataheader[f'{col}o{o}'] = self.dataord.loc[o, col]



    # def __repr__(self):
    #     return


    def __str__(self):
        # vars(self)
        return f''
    

    def plot_spectrum(self, ax=None, ords=None, x='w', y='f', wmin=None, wmax=None, normflux=None, offset=0, legend=False, legendloc=None, xunit='A', xlabel=None, ylabel='Flux', title='', lw=1, linestyle='-', alpha=1, alphaother=0.7, zorder=0, color=None, colorother=None, cmap=None, cbar=False):
        """Plot spectrum flux vs wavelength (or pixel), for the orders in `ords`.
        
        Parameters
        ----------
        x, y : str
            Parameters to plot from `Spectrum.dataspec`. Default to 'w' and 'f'.
        ords: list-like, optional
            Orders to plot. If None, plot all orders.
        wmin, wmax: floats, optional
            Wavelength range to plot. Must be in the same units as `w` (and `xunit`).

        normflux : function
            Function to normalise the spectrum. Default is None. max, mean...
        offset : float
            Offset to add to the spectrum. To ne used with `normflux` and pixel in the x axis. Default is None.
        c : str, optional
            Plot all orders with the same color, alternating alpha. Overrides `cmap`.
        cmap : str, optional
        cbar : bool, optional
        TODO add colorbar if using cmap
        """
        if ax is None: ax = plt.gca()
        if ords is None: ords = self.ords
        if np.issubdtype(type(ords), np.integer): ords = [ords]  # make sure it is a list

        # Wavelength range
        if wmin is None: wmin = np.nanmin(self.dataspec[x][ords])
        if wmax is None: wmax = np.nanmax(self.dataspec[x][ords])
        mp = (self.dataspec[x] >= wmin) & (self.dataspec[x] <= wmax)

        # Normalisation and offset
        if normflux is not None:
            yp = np.array([self.dataspec[y][o] / normflux(self.dataspec[y][o]) + offset*i for i, o in enumerate(self.ords)])
            print(yp.shape, self.dataspec[y].shape)
        else:
            yp = self.dataspec[y]

        # Colors
        if (cmap is not None) and (color is None):  # use cmap
            cmap = mpl.colormaps.get_cmap(cmap)
            norm = mpl.colors.Normalize(vmin=ords[0], vmax=ords[-1])
        elif (cmap is None) and (color is None):  # use default color cycle
            lencol = len(plt.rcParams['axes.prop_cycle'].by_key()['color'])
        
        a = alpha
        # Plot
        for o in ords:
            mpo = mp[o]
            if not any(mpo): continue  # skip if all pixels masked

            # Color and alpha
            if (cmap is not None) and (color is None):
                c = cmap(norm(o))
            elif (cmap is None) and (color is None):
                c = plt.rcParams['axes.prop_cycle'].by_key()['color'][o%lencol]
            elif (color is not None) and (colorother is None):
                c = color
                a = alpha if o % 2 == 0 else alphaother
            elif (color is not None) and (colorother is not None):
                c = color if o % 2 == 0 else colorother
                a = alpha if o % 2 == 0 else alphaother
            
            ax.plot(self.dataspec[x][o][mpo], yp[o][mpo], c=c, alpha=a, lw=lw, linestyle=linestyle, zorder=zorder, label=f'{o}')
        if legend: ax.legend(loc=legendloc)
        if xlabel is None: xlabel = wavelength_label(x=xunit)  # Assumes plotting wavelength in x-axis
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
    

    def plot_spectrum_pix(self, ax=None, ords=None, x='pix', y='f', pmin=None, pmax=None, normflux=np.nanmax, offset=1, legend=False, legendloc=None, xunit=None, xlabel='Pixel', ylabel='Flux', title='', lw=1, linestyle='-', alpha=1, zorder=0, color=None, cmap=None, cbar=False):
        """Plot spectrum normalised flux + offset vs pixel for the orders in `ords`.
        
        Parameters
        ----------
        x, y : str
            Parameters to plot from `Spectrum.dataspec`. Default to 'w' and 'f'.
        ords: list-like, optional
            Orders to plot. If None, plot all orders.
        wmin, wmax: floats, optional
            Wavelength range to plot. Must be in the same units as `w` (and `xunit`).

        normflux : function
            Function to normalise the spectrum. Default is None. max, mean...
        offset : float
            Offset to add to the spectrum. To ne used with `normflux` and pixel in the x axis. Default is None.
        c : str, optional
            Plot all orders with the same color, alternating alpha. Overrides `cmap`.
        cmap : str, optional
        cbar : bool, optional
        TODO add colorbar if using cmap
        """
        if ax is None: ax = plt.gca()
        if ords is None: ords = self.ords
        if np.issubdtype(type(ords), np.integer): ords = [ords]  # make sure it is a list

        # Wavelength range
        if pmin is None: pmin = np.nanmin(self.dataspec[x][ords])
        if pmax is None: pmax = np.nanmax(self.dataspec[x][ords])
        mp = (self.dataspec[x] >= pmin) & (self.dataspec[x] <= pmax)

        # Normalisation and offset
        if normflux is not None:
            yp = np.array([self.dataspec[y][o] / normflux(self.dataspec[y][o]) + offset*i for i, o in enumerate(self.ords)])
            print(yp.shape, self.dataspec[y].shape)
        else:
            yp = self.dataspec[y]

        # Colors
        if (cmap is not None) and (color is None):  # use cmap
            cmap = mpl.colormaps.get_cmap(cmap)
            norm = mpl.colors.Normalize(vmin=ords[0], vmax=ords[-1])
        elif (cmap is None) and (color is None):  # use default color cycle
            lencol = len(plt.rcParams['axes.prop_cycle'].by_key()['color'])
        
        a = alpha
        # Plot
        for o in ords:
            mpo = mp[o]
            if not any(mpo): continue  # skip if all pixels masked

            # Color and alpha
            if (cmap is not None) and (color is None):
                c = cmap(norm(o))
            elif (cmap is None) and (color is None):
                c = plt.rcParams['axes.prop_cycle'].by_key()['color'][o%lencol]
            ax.plot(self.dataspec[x][o][mpo], yp[o][mpo], c=c, alpha=a, lw=lw, linestyle=linestyle, zorder=zorder, label=f'{o}')
        if legend: ax.legend(loc=legendloc)
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title, xlim=(pmin, pmax))
        ax.minorticks_on()
        return ax
    

    def fig_spectrum_pix(self, filout=None, sh=False, sv=True, figsize=(16, 16), **kwargs):
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        ax = self.plot_spectrum_pix(ax=ax, **kwargs)
        if sh: plt.show()
        if sv:
            if filout is None: filout = f'{self.dirout}/{self.filname}_spec_pix.pdf'
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
            mpixs = np.ones_like(self.dataspec['f'][0], dtype=bool)
        else:
            mpixs = np.zeros_like(self.dataspec['f'][0], dtype=bool)
            mpixs[pixs] = True
        extent = [pixs[0], pixs[-1], ords[0], ords[-1]]
        # Flux limits
        if (vmin is not None) and (vmax is not None): extend = 'both'
        elif (vmin is not None) and (vmax is None): extend = 'min'
        elif (vmin is None) and (vmax is not None): extend = 'max'
        else: extend = None
        im = ax.imshow(self.dataspec['f'][mords,:][:,mpixs], origin='lower', interpolation='none', aspect='auto', extent=extent, vmin=vmin, vmax=vmax, cmap=cmap, rasterized=True)
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


    def plot_dataord(self, y, ax=None, x=None, z=None, ords=None, xlabel=None, ylabel=None, zlabel=None, title='', zorder=0, s=100, edgecolors='k', linewidths=0.5, ecolor=None, elinewidth=1.5, capsize=3, markeredgewidth=1.5, **kwargs):
        """
        kwargs for scatter plot (not errorbar)
        cmap in kwargs
        """
        if ax is None: ax = plt.gca()
        if x is None: x = self.ords
        if ords is None:
            m = np.ones_like(self.ords, dtype=bool)
        else:
            m = np.zeros_like(self.ords, dtype=bool)
            m[ords] = True
        zlabel = z if zlabel is None else zlabel
        dataz = self.dataord[z] if z is not None else None  # TODO Revise this, might fail if is None, or if z is not in dataord
        # Plot scatter
        sc = ax.scatter(x[m], self.dataord[y][m], c=dataz[m], zorder=0, s=s, edgecolors=edgecolors, linewidths=linewidths, **kwargs)
        # Error bars
        # self.dataord[y+'_err'] = np.ones(len(self.dataord[y])) * 5  # TESTING errobars
        if y+'_err' in self.dataord.columns:
            if ecolor is None: ecolor = mpl.colors.to_rgba('0.5', 0.8)
            ax.errorbar(x[m], self.dataord[y][m], yerr=self.dataord[y+'_err'][m], fmt='none', ecolor=ecolor, elinewidth=elinewidth, capsize=capsize, markeredgewidth=markeredgewidth, zorder=zorder-1)
        # Colorbar
        if z is not None:
            cbar = plt.colorbar(sc, ax=ax, label=zlabel)
            cbar.minorticks_on()
        # Style
        if xlabel is None: xlabel = 'Order'
        if ylabel is None: ylabel = y
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
        ax.minorticks_on()
        # Force x-labels to be integers
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        return ax


    def fig_dataord(self, y, filout=None, sh=False, sv=True, figsize=(8, 6), **kwargs):
        """
        """
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        # ax = self.plot_dataord(y, ax=ax, **kwargs)
        ax = self.plot_dataord(y, ax=ax, **kwargs)
        if sh: plt.show()
        if sv:
            if filout is None: filout = f'{self.dirout}/{self.filname}_ord_{y}.pdf'
            fig.savefig(filout)
        plt.close()
        return


class Spectra():
    """
    Collection of (time series) spectra of the same instrument.
    """

    def __init__(self, lisfil, inst, dirout='./', obj=None, tag='', deleteindividual=False, ext=0,
    # carmnir parameters
    ordcut=True, saveordnoncut=False,
    # harps and harpsn parameters
    readblaze=True, dirblaze=None, filblaze=None,
    ):
        """
        tag : str, optional
            Tag to add to the object name (e.g. type of reduction, sky subtracted...).

        ext : int, optional
            Header extension to read.

        Accessing the properties of Spectra
        ```
        # `lisfil` is a list of 4 spectra with 61 orders and 4096 pixels each
        
        # Acces the flux `f` of all spectra
        >>> spectra = Spectra(lisfil, inst)
        >>> spectra.dataspec['f'].shape
        (4, 61, 4096)
        
        # Access the flux `f` of order 4 for all spectra:
        >>> spectra.dataspec['f'][:,4,:].shape
        (4, 4096)

        # Access the flux `f` of orders 4 to 8 and pixels 1500 to 1600 for the first 3 observations:
        >>> spectra.dataspec['f'][:3,4:9,1500:1601].shape
        (3, 5, 101)
        ```
        """
        self.lisfil = lisfil
        self.lisfil_initial = lisfil
        self.inst = inst
        self.obj = obj
        self.tag = tag
        self.dirout = dirout

        if not os.path.exists(self.dirout): os.makedirs(self.dirout)

        headertable = read_header_kw_table()

        # Read files into Spectrum objects
        self.lisspec, self.lisfil = read_spectra(self.lisfil_initial, self.inst, returnclass=True, dirout=dirout, ordcut=ordcut, saveordnoncut=saveordnoncut, headertable=headertable)

        # Get filenames
        self.lisfilname = [spec.filname for spec in self.lisspec]

        # Rearrange individial data into dictionary of properties for all spectra in `dataspec`
        self.dataspec = {}
        for k in self.lisspec[0].dataspec.keys():
            self.dataspec[k] = np.array([spec.dataspec[k] for spec in self.lisspec])

        # For general properties shared by all spectra, add only the value from the first spectrum
        self.nord = self.lisspec[0].nord
        self.ords = self.lisspec[0].ords
        self.ords_real = self.lisspec[0].ords_real
        self.ord_ref = self.lisspec[0].ord_ref
        self.pixel_ms = self.lisspec[0].pixel_ms
        self.npix = self.lisspec[0].npix
        if self.obj is None: self.obj = self.lisspec[0].obj
        if self.tag is None: self.tag = self.lisspec[0].tag

        # Get parameters from header in `dataheader` attribute
        # Transform list of dictionaries into pandas dataframe
        lisdataheader = [sp.dataheader for sp in self.lisspec]
        self.dataheader = pd.DataFrame(lisdataheader, index=self.lisfilname)
        del lisdataheader

        # Get per order parameters in `dataord` attribute
        # dataord is a dictionary with a pandas dataframe for each propety, e.g. dataord['snr'] is a pandas dataframe with the S/N of all orders (columns) for all spectra (index). Can shift orders and spectra with `dataord['snr'].T`.
        lisdataord = [sp.dataord for sp in self.lisspec]
        self.dataord = {}
        for k in lisdataord[0].keys():
            self.dataord[k] = []
            for i, d in enumerate(lisdataord):
                self.dataord[k].append(d[k])
            self.dataord[k] = pd.DataFrame(self.dataord[k], index=self.lisfilname)
        
        # Add per order parameters from `dataord` to `dataheader`. Keywords will be '{parameter}o{ord}' where parameter is the column in `dataord` (e.g. 'snr') and ord is the order number.
        for k in self.dataord.keys():
            # Update colum names from order number to {k}o{ord}
            rename_cols = {o: f'{k}o{o}' for o in self.dataord[k].columns}
            dftemp = self.dataord[k].rename(columns=rename_cols)
            self.dataheader = pd.concat([self.dataheader, dftemp], axis=1)

        # Delete individual Spectrum objects to save memory
        if deleteindividual: del self.lisspec

    # def __repr__(self):
    #     return


    def __str__(self):
        # vars(self)
        return f''
    

    def plot_spectra(self, ax=None, ords=None, lisspec=None, wmin=None, wmax=None, legendlabel=None, legendwhich='first', legendloc=None, xunit='A', xlabel=None, ylabel='Flux', title='', cprop=None, cprop_all=False, cbar=None, cbarlabel=None, cmap=None, lw=1, linestyle0='-', linestyle1='-', alpha0=1, alpha1=0.7, zorder=1):
        """
        Plot spectra flux vs wavelength, orders in `ords`.

        Plot different observations with different colors, and optionally alternating alpha for different orders.

        Parameters
        ----------
        ords : list-like, optional
            List of orders to plot. If `None` (default), all orders are plotted.
        
        lisspec : list-like, optional
            List of observations (the number) to plot. Default is None, which plots all spectra in `lisspec`.

        cprop : array, optional
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
        
        legendwhich : 'first' or 'all', optional
            To label (with `legendlabel`) only the first spectrum in the list, or all spectra in the list. Default is 'first'.
        legendlabel : str, optional, default None
            Label for the legend. Can label the first order of a single spectrum in the list, or the first order of all spectra in the list.
            If only the first spectrum is labelled (with `legendwhich='first'`), the labe should be somethign like 'Observations' or 'Spectra'.
            If all spectra are labelled, the label will be the name of each spectra.
        """
        if ax is None: ax = plt.gca()
        if lisspec is None: lisspec = np.arange(len(self.lisspec))
        if np.issubdtype(type(lisspec), np.integer): lisspec = [lisspec]  # make sure it is a list
        if ords is None: ords = self.ords
        if np.issubdtype(type(ords), np.integer): ords = [ords]  # make sure it is a list
        
        # Wavelength range
        if wmin is None: wmin = np.nanmin(self.dataspec['w'][lisspec][:,ords])
        if wmax is None: wmax = np.nanmax(self.dataspec['w'][lisspec][:,ords])
        mp = (self.dataspec['w'] >= wmin) & (self.dataspec['w'] <= wmax)

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
        for i, idx in enumerate(lisspec):
            spec = self.lisspec[idx]
            if cprop is not None:
                c = cmap(norm(cprop[idx]))
            elif cmap is not None:
                c = cmap(norm(i))
            else:
                c = plt.rcParams['axes.prop_cycle'].by_key()['color'][idx]
            # c = cmap(norm(i)) if cmap is not None else plt.rcParams['axes.prop_cycle'].by_key()['color'][idx]
            # Loop over orders
            for o in ords:
                mpo = mp[idx][o]
                if not any(mpo): continue  # skip if all pixels masked
                a = alpha0 if o % 2 == 1 else alpha1
                linestyle = linestyle0 if o % 2 == 1 else linestyle1
                if legendwhich == 'all':
                    label = f'{spec.filname}' if (o == ords[0]) else None
                elif legendwhich == 'first':
                    label = legendlabel if (o == ords[0]) and (i == 0) else None
                ax.plot(spec.dataspec['w'][o][mpo], spec.dataspec['f'][o][mpo], c=c, alpha=a, lw=lw, linestyle=linestyle, label=label, zorder=zorder)
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
        TODO
        """
        pass
        return


    def plot_dataord(self, y, ax=None, ords=None, lisspec=None, z=None, xlabel='Order', ylabel=None, zlabel=None, title='', zorder=0, s=100, edgecolors='k', linewidths=0.5, ecolor=None, elinewidth=1.5, capsize=3, markeredgewidth=1.5, cmap=None, **kwargs):
        """
        Plot property `y` vs orders for orders in `ords` and spectra in `lisspec`. Optional color code by general spectrum property `z` from self.dataheader.

        Expect lisspec and ord to be sorted (increasing).
        """
        if ax is None: ax = plt.gca()
        if lisspec is None: lisspec = np.arange(len(self.lisfil))
        if np.issubdtype(type(lisspec), np.integer): lisspec = [lisspec]  # make sure it is a list
        if ords is None: ords = self.ords
        if np.issubdtype(type(ords), np.integer): ords = [ords]  # make sure it is a list

        # if x is None: x = self.ords
        # if ords is None:
        #     m = np.ones_like(self.ords, dtype=bool)
        # else:
        #     m = np.zeros_like(self.ords, dtype=bool)
        #     m[ords] = True

        # Get data from dataframe and cut by spectra (y and z) and orders (y)
        datay = self.dataord[y].iloc[lisspec, ords]
        # self.dataord[y+'_err'] = self.dataord[y].copy() / 5  # TESTING errobars
        if y+'_err' in self.dataord.keys(): datayerr = self.dataord[y+'_err'].iloc[lisspec, ords]


        zlabel = z if zlabel is None else zlabel
        dataz = self.dataheader[z].iloc[lisspec] if z is not None else None

        # Colors
        if (z is None) and (cmap is None):  # use default color cycle
            pass
        elif (z is None) and (cmap is not None):  # use cmap with spectrum index
            cmap = mpl.colormaps.get_cmap(cmap)
            norm = mpl.colors.Normalize(vmin=np.nanmin(lisspec), vmax=np.nanmax(lisspec))
        elif (z is not None):  # use cmap with property z
            if cmap is None: cmap = 'viridis'  # default cmap
            cmap = mpl.colormaps.get_cmap(cmap)
            norm = mpl.colors.Normalize(vmin=np.nanmin(dataz), vmax=np.nanmax(dataz))

        # Plot scatter
        for i, _ in enumerate(datay.index):
            if z is not None:
                c = cmap(norm(dataz.iloc[i]))
            elif cmap is not None: c = cmap(norm(i))
            else: c = None
            sc = ax.scatter(ords, datay.iloc[i], c=c, zorder=zorder, s=s, edgecolors=edgecolors, linewidths=linewidths, **kwargs)
            # Error bars
            if y+'_err' in self.dataord.keys():
                if ecolor is None: ecolor = mpl.colors.to_rgba('0.5', 0.8)
                ax.errorbar(ords, datay.iloc[i], yerr=datayerr.iloc[i], fmt='none', ecolor=ecolor, elinewidth=elinewidth, capsize=capsize, markeredgewidth=markeredgewidth, zorder=zorder-1)
        # Colorbar
        if z is not None:
            cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label=zlabel)
            cb.minorticks_on()
        # Style
        if ylabel is None: ylabel = y
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
        ax.minorticks_on()
        # Force x-labels to be integers
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        return ax


    def fig_dataord(self, y, filout=None, sh=False, sv=True, figsize=(8, 5), **kwargs):
        """
        """
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        ax = self.plot_dataord(y, ax=ax, **kwargs)
        if sh: plt.show()
        if sv:
            if filout is None: filout = f'{self.dirout}/{self.filname}_ord_{y}.pdf'
            fig.savefig(filout)
        plt.close()
        return





