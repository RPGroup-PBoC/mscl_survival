import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import skimage.morphology
import skimage.measure


def contour_seg(image, level, selem='default'):
    """
    Identifies contours around dark objects in a phase contrast image.

    Parameters
    ----------
    image: 2d-array
        Phase contrast image of interest. This shoul
    level: float
        Level at which to draw contours on black top-hat filtered image.
        Default value is 4000. This is dependent on the image, but should
        be robust for images of similar illumination.
    selem: 2d-array or string
        Structuring element to use for the black top-hat filtering procedure
        Default value is a disk with a diameter of 20 pixels.

    Returns
    -------
    conts : 1d-array
        List of contour coordinates. Each entry of this array comes as
        an x,y pair of arrays. Has the same length as the number of
        contoured objects

    """

    # Apply the white top-hat filter.
    if selem == 'default':
        selem = skimage.morphology.disk(20)
    im_filt = skimage.morphology.black_tophat(image, selem)

    # Find the contours and return.
    conts = skimage.measure.find_contours(im_filt, level)
    return conts


def set_plotting_style(return_colors=True):
    """
    Sets the plotting style.

    Parameters
    ----------
    return_colors: Bool
        If True, this will also return a palette of eight color-blind safe
        colors with the hideous yellow replaced by 'dusty purple.'
    """
    rc = {'lines.linewidth': 2,
          'axes.facecolor': '#E3DCD0',
          'xtick.labelsize': 'large',
          'font.family': 'Lucida Sans Unicode',
          'grid.linestyle': ':',
          'grid.linewidth': 0.85,
          'grid.color': '#ffffff',
          'mathtext.fontset': 'stixsans',
          'mathtext.sf': 'sans',
          'legend.frameon': True}
    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    sns.set_style('darkgrid', rc=rc)
    sns.set_context('paper', rc=rc)
    colors = sns.color_palette('colorblind', n_colors=8)
    colors[4] = sns.xkcd_palette(['dusty purple'])[0]

    if return_colors:
        return colors
