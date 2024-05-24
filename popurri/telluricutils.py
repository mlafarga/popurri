import os
from astropy.io import fits
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

from . import peakutils
from . import plotutils
from .plotutils import wavelength_label

# Constants
C_MS = 2.99792458 * 1.e8  # Light speed [m/s]
C_KMS = 2.99792458 * 1.e5  # Light speed [km/s]

# Strings
angstromstr = r'$\mathrm{\AA}$'
angstromstrunit = r'$[\mathrm{\AA}]$'


###############################################################################

# Telluric models
# ---------------

# ESO SkyCalc model telluric spectrum
# -----------------------------------

def read_esoskycalc_telluric(filin, extra=False):
    """
    Read transmision telluric model spectrum generated with ESO SkyCalc.

    https://www.eso.org/observing/etc/skycalc/

    Notes
    -----
    
    FITS file structure:
        >>> hdu.info()
        Filename: ./esoskycalc/lasilla_w300-1000_R100000_FWHM3.20bin_vac.fits
        No.    Name      Ver    Type      Cards   Dimensions   Format
          0  PRIMARY       1 PrimaryHDU      67   ()
          1                1 BinTableHDU    104   120398R x 18C   ['1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D']

    FITS file data from extension 0 header:
        COMMENT column lam:     Vacuum wavelength in nm
        COMMENT column flux:    sky emission radiance flux in ph/s/m2/micron/arcsec2
        COMMENT column trans:   sky transmission

        COMMENT Individual emission components in ph/s/m2/micron/arcsec2
        COMMENT column flux_sml:    scattered moonlight
        COMMENT column flux_ssl:    scattered starlight
        COMMENT column flux_zl:     zodiacal light
        COMMENT column flux_tie:    telescope/instrument thermal emission
        COMMENT column flux_tme:    molecular emission lower atmosphere
        COMMENT column flux_ael:    airglow emission Lines
        COMMENT column flux_arc:    airglow/residual emission

        COMMENT Individual transmission components
        COMMENT column trans_o3:    ozone uv/optical absorption
        COMMENT column trans_ma:    molecular absorption
        COMMENT column trans_rs:    rayleigh scattering
        COMMENT column trans_ms:    mie scattering
    """
    data = {}
    with fits.open(filin) as hdu:
        data['w'] = hdu[1].data['lam']  # Vacuum wavelength in nm
        data['w'] = data['w'] * 10  # nm -> A
        data['femission'] = hdu[1].data['flux']  # sky emission radiance flux in ph/s/m2/micron/arcsec2
        data['ftransmission'] = hdu[1].data['trans']  # sky transmission
        if extra == True:
            dataextra = {
                'trans_o3': hdu[1].data['trans_o3'], 'trans_ma': hdu[1].data['trans_ma'], 'trans_rs': hdu[1].data['trans_rs'], 'trans_ms': hdu[1].data['trans_ms'],
                'flux_sml': hdu[1].data['flux_sml'], 'flux_ssl': hdu[1].data['flux_ssl'], 'flux_zl': hdu[1].data['flux_zl'], 'flux_tie': hdu[1].data['flux_tie'], 'flux_tme': hdu[1].data['flux_tme'], 'flux_ael': hdu[1].data['flux_ael'], 'flux_arc': hdu[1].data['flux_arc']
            }
        else:
            dataextra = None
        data.update(dataextra)
    # return w, ftransmission, femission, lisextra
    return data


def read_esoskycalc_telluric_transmission(filin, extra=False):
    """OLD FUNCTION!

    Read transmision telluric model spectrum generated with ESO SkyCalc.

    https://www.eso.org/observing/etc/skycalc/

    Notes
    -----
    
    FITS file structure:
        >>> hdu.info()
        Filename: ./esoskycalc/lasilla_w300-1000_R100000_FWHM3.20bin_vac.fits
        No.    Name      Ver    Type      Cards   Dimensions   Format
          0  PRIMARY       1 PrimaryHDU      67   ()
          1                1 BinTableHDU    104   120398R x 18C   ['1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D', '1D']

    FITS file data from extension 0 header:
        COMMENT column lam:     Vacuum wavelength in nm
        COMMENT column flux:    sky emission radiance flux in ph/s/m2/micron/arcsec2
        COMMENT column trans:   sky transmission

        COMMENT Individual emission components in ph/s/m2/micron/arcsec2
        COMMENT column flux_sml:    scattered moonlight
        COMMENT column flux_ssl:    scattered starlight
        COMMENT column flux_zl:     zodiacal light
        COMMENT column flux_tie:    telescope/instrument thermal emission
        COMMENT column flux_tme:    molecular emission lower atmosphere
        COMMENT column flux_ael:    airglow emission Lines
        COMMENT column flux_arc:    airglow/residual emission

        COMMENT Individual transmission components
        COMMENT column trans_o3:    ozone uv/optical absorption
        COMMENT column trans_ma:    molecular absorption
        COMMENT column trans_rs:    rayleigh scattering
        COMMENT column trans_ms:    mie scattering
    """
    with fits.open(filin) as hdu:
        w = hdu[1].data['lam']  # Vacuum wavelength in nm
        w = w * 10  # nm -> A
        # femission = hdu[1].data['flux']  # sky emission radiance flux in ph/s/m2/micron/arcsec2
        ftransmission = hdu[1].data['trans']  # sky transmission
        if extra == True:
            lisextra = {'trans_o3': hdu[1].data['trans_o3'], 'trans_ma': hdu[1].data['trans_ma'], 'trans_rs': hdu[1].data['trans_rs'], 'trans_ms': hdu[1].data['trans_ms']}
        else:
            lisextra = None
    return w, ftransmission, lisextra


# Model class
# -----------

class ModelEsoSkyCalc():
    """
    w : array
        Wavelength [A].
    ftransmission, femission : array
        Transmission and emission flux.
    extra: Dictionary with keys 'trans_o3', 'trans_ma','trans_rs', 'trans_ms'
    """
    def __init__(self, filin, dirout='./', diroutplot='./', tag='', wk='w', fk='ftransmission'):
        self.filin = filin
        self.tag = tag
        self.utag = '_' + tag  # tag with underscore for file names
        self.dirout = dirout
        self.diroutplot = diroutplot
        for dout in [self.dirout, self.diroutplot]:
            if not os.path.exists(dout): os.makedirs(dout)

        self.wk = wk
        self.fk = fk

        # Read model from filin
        self.data = read_esoskycalc_telluric(self.filin, extra=True)
        self.wunit = 'A'

        # Attributes to be set by running below functions
        self.imin = None
        self.imaxl = None
        self.imaxr = None
        self.npixmin = None
        self.maskpix = None
        self.depthmin = None
        self.maskshallow = None
        self.fcut = None
        self.dw = None
        self.masks = {}


    def cut_w(self, wmin, wmax):
        """Cut model to a specific wavelength range."""
        w = self.data[self.wk]
        # Adjust wavelength (good = True in mask)
        if (wmin is not None) and (wmax is not None):
            mask = (w >= wmin) & (w <= wmax)
        elif (wmin is not None) and (wmax is None):
            mask = w >= wmin
        elif (wmin is None) and (wmax is not None):
            mask = w <= wmax
        elif (wmin is None) and (wmax is None):
            mask = w >= 0.
        for k in self.data.keys():
            self.data[k] = self.data[k][mask]


    def find_lines(self):
        """Find minima in the telluric spectrum."""
        self.imin = peakutils.idxmin_custom(self.data[self.fk])
        self.imaxl, self.imaxr = peakutils.idxminends_custom(self.data[self.fk], self.imin)
        return


    def rmv_lines_npix(self, npixmin=3):
        """Find and remove lines with less than npixmin pixels."""
        self.npixmin = npixmin

        # Make sure find_lines has been run, or if not, run it
        if self.imin is None:
            self.find_lines()

        maskpix = []  # lines to keep
        for i, imin_i in enumerate(self.imin):
            if self.imaxr[i] - self.imaxl[i] > npixmin:
                maskpix.append(i)
        self.maskpix = maskpix

        # Update
        self.imin, self.imaxl, self.imaxr = self.imin[maskpix], self.imaxl[maskpix], self.imaxr[maskpix]
        return


    def rmv_lines_depth(self, depthmin=1-5*1e-5):
        """Find and remove lines with depth less than depthmin (i.e. shallow lines)."""
        self.depthmin = depthmin

        # Make sure find_lines has been run, or if not, run it
        if self.imin is None:
            self.find_lines()
        maskshallow = self.data[self.fk][self.imin] <= depthmin
        self.maskshallow = maskshallow

        # Update
        self.imin, self.imaxl, self.imaxr = self.imin[maskshallow], self.imaxl[maskshallow], self.imaxr[maskshallow]
        return


    # --- Other functions ---
    # Should test in a notebook
    #
    # # Seems to work as well as the custom method (i.e. gets the same minima)
    # # Need to separate imax into left and right
    # imin, minprops = find_peaks(1 - f, distance=1)
    # imax, maxprops = find_peaks(f, distance=1)
    # # Can select peaks within the function already
    # imin1, minprops = find_peaks(1 - f, distance=1, threshold=10**(-6), width=(None, None), plateau_size=(None, None))
    # imax1, maxprops = find_peaks(f, distance=1, threshold=10**(-6), width=(None, None), plateau_size=(None, None))
    #
    # # argelmin / argrelmax: Doesn't work well with flat minima/maxima
    # imin = argrelmin(f)[0]
    # imax = argrelmax(f)[0]

    # --- Different approach ---
    # ------------------
    # # ISSUE with this method: do not know the amplitude of the peaks, so cannot make a proper mask
    # # Select telluric lines (good = True in mask)
    # mask = f >= args.fcut
    # # Make flux good = 0 and bad = 1
    # fmask = np.ones_like(f)
    # fmask[mask] = 0
    #
    # # Plot
    # fig, ax = plt.subplots(figsize=(16, 8))
    # ax.plot(w, f, marker='.', rasterized=True, alpha=0.8)
    # ax.plot(w[mask], f[mask], marker='.', rasterized=True, alpha=0.8)
    # ax = telluricutils.plot_mask(w, fmask, ax=ax)
    # ax.set_ylabel('Flux')
    # ax.set_xlabel('Wavelength ' + angstromstrunit)
    # ax.minorticks_on()
    # plt.tight_layout()
    # plt.show()
    # plt.close()
    # ------------------

    def make_mask(self, fcut=0.999, dw=0.001, sv=True, filout=None, verbose=True, svobj=False):
        """
        Make telluric mask.
        
        svojb : bool
            Save mask in object dictionary (as a Mask object), with fcut as key. Useful if making several masks from the same model and want to compare them.
        """
        self.fcut = fcut
        self.dw = dw
        w =  self.data[self.wk]
        f =  self.data[self.fk]
        imin = self.imin
        imaxl = self.imaxl
        imaxr = self.imaxr

        if imin is None:
            self.find_lines()
            imin = self.imin
            imaxl = self.imaxl
            imaxr = self.imaxr

        # Mask shallow lines
        m = f[imin] <= fcut
        imin_cut = imin[m]
        imaxl_cut = imaxl[m]
        imaxr_cut = imaxr[m]

        # Merge overlapping lines
        wt1, wt2 = join_overlap_wlimits(w[imaxl_cut], w[imaxr_cut])

        # Make telluric mask (0110)
        #   I.e. each line is:
        #       w: imaxl_i-dw imaxl_i imaxr_i imaxr_i+dw
        #       f: 0          1       1       0
        wt, ft = wlimits2mask(wt1, wt2, dw=dw)

        # Output
        if sv:
            textextra = f'_npixmin{self.npixmin}' if self.npixmin is not None else ''
            if filout is None: filout = os.path.join(self.dirout, f'mask{self.utag}_fcut{fcut:.6f}{textextra}.mas')
            np.savetxt(filout, np.vstack([wt, ft]).T, fmt='%.8f')
            if verbose: print('Telluric mask saved:', filout)
        
        # Save object
        # if svobj: self.masks[fcut] = {'w': wt, 'f': ft, 'w1': wt1, 'w2': wt2}
        if svobj: 
            self.masks[fcut] = Mask(filin=None, w=wt, f=ft, dirout=self.dirout, diroutplot=self.diroutplot, tag=self.tag, join_overlap=False, broaden_dv=None)
        return
    
    # =========================================================================

    # Plots

    def plot_spec(self, ax=None, xunit='A', color='k', xlabel=None, ylabel='Flux', leglabel='', title='', zorder=100, **kwargs):
        if ax is None: ax = plt.gca()
        # ax.plot(self.data[self.wk], self.data[self.fk], color=color, label=leglabel, **kwargs)
        ax.plot(self.data[self.wk], self.data[self.fk], color=color, zorder=zorder, **kwargs)
        if xlabel is None: xlabel = wavelength_label(x=xunit)
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
        return ax
    

    def fig_spec(self, lisfk=None, filout=None, sh=False, sv=True, figsize=(16, 4), svext=['pdf'], **kwargs):
        if lisfk is None:
            lisfk = self.data.keys()
            lisfk.remove(self.wk)
        if isinstance(lisfk, str): lisfk = [lisfk]
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        for fk in lisfk:
            ax = self.plot_spec(ax=ax, leglabel=fk, **kwargs)
        if filout is None:
            filout = os.path.join(self.diroutplot, 'spec' + self.utag)
        plotutils.figout(fig, sv=sv, filout=filout, svext=svext, sh=sh)
        return


    def plot_lines(self, ax=None, ylim=None, fimin={'color':'0.5', 'zorder':-10}, fimaxl={'color':'b', 'linestyle':'dashed', 'alpha':0.5, 'zorder':-10}, fimaxr={'color':'r', 'linestyle':'dotted', 'alpha':0.5, 'zorder':-10}):
        """
        Plot identified lines.
        
        It is recommended to do `plot_spec()` first, get ylimits of tge plots and pass them here with `ylim`.
        """
        if ax is None: ax = plt.gca()
        w = self.data[self.wk]
        if ylim is not None: y0, y1 = ylim[0], ylim[1]
        else: y0, y1 = 0, 1
        # Plot lines
        ax.vlines(w[self.imin], y0, y1, **fimin)
        ax.vlines(w[self.imaxl], y0, y1, **fimaxl)
        ax.vlines(w[self.imaxr], y0, y1, **fimaxr)
        # ax.set_title(f'Removed lines with <= {self.npixmin} pixels and depth >= {self.depthmin:.6f}')
        return ax


    def fig_spec_lines(self, filout=None, sh=False, sv=True, figsize=(16, 4), svext=['pdf'], fimin={'color':'0.5', 'zorder':-10}, fimaxl={'color':'b', 'linestyle':'dashed', 'alpha':0.5, 'zorder':-10}, fimaxr={'color':'r', 'linestyle':'dotted', 'alpha':0.5, 'zorder':-10}, **kwargs):
        """
        kwargs : dict
            Only passed to `plot_spec()`.
        """
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        # Spec
        ax = self.plot_spec(ax=ax, leglabel=self.fk, **kwargs)
        # Lines
        ylim = ax.get_ylim()
        ax = self.plot_lines(ax=ax, ylim=ylim, fimin=fimin, fimaxl=fimaxl, fimaxr=fimaxr)
        ax.set_ylim(ylim)
        # Format
        ax.minorticks_on()
        # Save
        if filout is None:
            textextra0 = f'_npixmin{self.npixmin}' if self.npixmin is not None else ''
            textextra1 = f'_depthmin{self.depthmin:.6f}' if self.depthmin is not None else ''
            textextra = textextra0 + textextra1
            filout = os.path.join(self.diroutplot, f'spec_lines{self.utag}{textextra}')
        plotutils.figout(fig, sv=sv, filout=filout, svext=svext, sh=sh)
        return


###############################################################################

# Telluric mask utils
# -------------------

def read_mask(filin, xunit='A'):
    """Read telluric mask from file.
    
    Columns: 0) wavelength, 1) flux.
    Telluric lines are indicated by 4 lines with flux "0110".
    Assumes that wavelength is in Angstroms.
    """
    w, f = np.loadtxt(filin, usecols=[0, 1], unpack=True)
    return w, f


def mask2wlimits(w, f):
    """Get wavelength limits of the telluric lines masked.

    If a line is given by:
        wa 0
        wb 1
        wc 1
        wd 0
    then get wb as the first limit (store it in `w1`) and wc as the second limit (store it in `w2`).
    """
    w1, w2 = [], []
    for i in range(len(w)-1):
        # Normal line 0110
        if f[i-1] == 0 and f[i] == 1 and f[i+1] == 1 and f[i+2] == 0:
            w1.append(w[i])
            w2.append(w[i+1])
        # Blended lines 01110
        elif f[i-1] == 0 and f[i] == 1 and f[i+1] == 1 and f[i+2] == 1 and f[i+3] == 0:
            w1.append(w[i])
            w2.append(w[i+2])
        # Blended lines 011110
        elif f[i-1] == 0 and f[i] == 1 and f[i+1] == 1 and f[i+2] == 1 and f[i+3] == 1 and f[i+4] == 0:
            w1.append(w[i])
            w2.append(w[i+3])
    w1, w2 = np.asarray(w1), np.asarray(w2)
    return w1, w2


def wlimits2mask(w1, w2, dw=0.001):
    """Convert telluric regions indicated by wavelengths w1 and w2 to a mask with format: flux = 0110 per region."""
    w, f = [], []
    for i in range(len(w1)):  # for each line
        w += [w1[i]-dw, w1[i], w2[i], w2[i]+dw]
        f += [0., 1., 1., 0.]
    w, f = np.asarray(w), np.asarray(f)
    return w, f


def is_there_overlap_wlimits(w1, w2):
    """Check if there is overlap between 2 consecutive telluric regions.
    
    Uses w limits w1 and w2 as inputs, not the raw mask data w and f.
    """
    ret = False
    for i in range(len(w1)-1):
        if w2[i] >= w1[i+1]:
            ret = True
            break
    return ret


def join_overlap_wlimits_once(w1, w2):
    """Join consecutive telluric regions if there is overlap.
    Only check once.

    Uses w limits w1 and w2 as inputs, not the raw mask data w and f.
    """
    w1new, w2new = [], []
    iused = -1
    for i in range(len(w1)-1):  # for each telluric region
        if i == iused:
            continue
        w1new.append(w1[i])
        if w2[i] >= w1[i+1]:  # Join lines
            # w1new.append(w1[i])
            w2new.append(w2[i+1])
            iused = i+1
        else:
            w2new.append(w2[i])

    return w1new, w2new


def join_overlap_wlimits(w1, w2):
    """Join consecutive telluric regions if there is overlap.
    Join until there is no overlap.
    """
    while is_there_overlap_wlimits(w1, w2):
        w1, w2 = join_overlap_wlimits_once(w1, w2)
    return w1, w2


def interp_mask_inverse(w, f, kind='linear'):
    """Interpolation function of the telluric mask inverted (so that lines are zeros, i.e. False) to be used to mask telluric regions in the data.

    Example
    -------

    >>> # wt, ft : Telluric mask
    >>> # ws, fs : Spectrum
    >>> # Make the inverted mask function
    >>> MaskInv, _ = interp_mask_inverse(wt, ft) # function
    >>> # Make the spectrum flux 0 where there are tellurics
    >>> fs_clean = MaskInv(ws)
    >>> # Plot
    >>> plt.plot(ws, fs, label='spec')
    >>> plt.plot(ws, fs_clean, label='spec clean')
    >>> plt.fill_between(wt, 0, ft, label='tell mask')
    """
    f_inv = np.array(~np.array(f, dtype=bool), dtype=float)
    MaskInv = interp1d(w, f_inv, kind=kind)
    return MaskInv, f_inv


def broaden_mask(w, f, dv):
    """Broaden mask by a velocity `dv` [m/s].
    
    w_broaden = w*(1+dv/(2.99792458*10**8))

    Uses mask format w and f, where one lines has flux 0110.

    Parameters
    ----------
    w, f : array
    dv : float
        Velocity to add/subtract to each telluric line, in  m/s.

    Returns
    -------
    w_broaden : array
    """

    C_MS = 2.99792458e8  # [m/s]

    dv = abs(dv)

    w_broaden = [[]]*len(w)
    idone = np.zeros(len(w))
    i = 0
    while i < len(w):
        if idone[i] == 0 and i < len(w)-1:
            # Line (where flux = 1)
            if f[i] == 1 and f[i+1] == 1:
                w_broaden[i] = w[i]*(1-dv/C_MS)
                w_broaden[i+1] = w[i+1]*(1+dv/C_MS)
                idone[i], idone[i+1] = 1, 1
                i = i+2
            # Base of the lines (where flux = 0)
            elif f[i] == 0 and f[i+1] == 0:
                w_broaden[i] = w[i]*(1+dv/C_MS)
                w_broaden[i+1] = w[i+1]*(1-dv/C_MS)
                idone[i], idone[i+1] = 1, 1
                i = i+2
            # Half base (can only happend at the beginning)
            elif f[i] == 0 and f[i+1] == 1:
                w_broaden[i] = w[i]
                i = i+1

        # Last value
        elif idone[i] == 0 and i == len(w)-1:
            if f[i] == 1:  # ...01
                w_broaden[i] = w[i]*(1-dv/C_MS)
            if f[i] == 0:  # ...10
                w_broaden[i] = w[i]*(1+dv/C_MS)
            i = i+1
    w_broaden = np.asarray(w_broaden)
    return w_broaden


def broaden_wlimits(w1, w2, dv):
    """Same as `broaden_mask`, but takes as input the wavelength limits of each mask line, instead of the wavelength and flux (0110).
    """
    C_MS = 2.99792458e8  # [m/s]

    dv = abs(dv)
    w1_broaden = np.asarray(w1)*(1-dv/C_MS)
    w2_broaden = np.asarray(w2)*(1+dv/C_MS)
    return w1_broaden, w2_broaden


# Plots


###############################################################################

# Classes
# -------

class Mask():
    """
    Telluric mask class.

    Can create from data in file `filin`, or from data already loaded with `w` and `f`. If both `filin` and `w` and `f` are given, `filin` has preference and is the on that will be used.
    """
    def __init__(self, filin=None, w=None, f=None, dirout='./', diroutplot='./', tag='', join_overlap=True, broaden_dv=None):
        self.filin = filin
        self.dirout = dirout
        self.diroutplot = diroutplot
        for dout in [self.dirout, self.diroutplot]:
            if not os.path.exists(dout): os.makedirs(dout)
        self.tag = tag
        self.utag = '_' + tag  # tag with underscore for file names
        self.join_overlap = join_overlap
        self.broaden_dv = broaden_dv  # [m/s]

        # Read mask or get it from w and f
        if self.filin is not None:
            self.w, self.f = read_mask(self.filin)
        else:
            if (w is None) or (f is None):
                raise ValueError('If filin is None, w and f must be given.')
            self.w, self.f = w, f

        # Broaden telluric lines
        if self.broaden_dv is not None:
            self.w = broaden_mask(self.w, self.f, self.broaden_dv)

        # Get wavelength limits of the lines in the mask
        self.w1, self.w2 = mask2wlimits(self.w, self.f)

        # Join lines that overlap
        if self.join_overlap:
            self.w1, self.w2 = join_overlap_wlimits(self.w1, self.w2)

        # Make mask to be used with data
        self.Mask, _ = interp_mask_inverse(self.w, self.f, kind='linear')


    def plot_mask(self, ax=None, fscale=1, zorder=-1, xunit='A', xlabel=None, ylabel='Flux', leglabel='', title='', color='0.5', alpha=0.6, **kwargs):
        if ax is None: ax = plt.gca()
        # if fscale != 1:
        #     f = f * fscale
        # ax.fill_between(w, f, zorder=zorder, label=leglab, **kwargs)
        ax.fill_between(self.w, self.f*fscale, zorder=zorder, color=color, alpha=alpha, label=leglabel, **kwargs)
        if xlabel is None: xlabel = wavelength_label(x=xunit)  # Assumes plotting wavelength in x-axis
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
        return ax


    def fig_mask(self, filout=None, sh=False, sv=True, svext=['pdf'], figsize=(16, 4), **kwargs):
        fig, ax = plt.subplots(figsize=figsize)
        ax = plot_mask(ax=ax, **kwargs)
        if filout is None:
            filout = os.path.join(self.diroutplot, 'mask' + self.utag)
        plotutils.figout(fig, sv=sv, filout=filout, svext=svext, sh=sh)
        return



