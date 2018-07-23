import os
import skimage.io
import skimage.measure
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import bokeh.plotting
import bokeh.io
import bokeh.palettes
import bokeh.layouts


# ---------------------------------------------------------------------------
# Plotting styles
# ---------------------------------------------------------------------------
def set_plotting_style(return_colors=True):
    """
    Sets the plotting style.

    Parameters
    ----------
    return_colors: Bool
        If True, this will also return a palette of eight color-blind safe
        colors with the hideous yellow replaced by 'dusty purple.'
    """
    rc = {'axes.facecolor': '#E3DCD0',
          'font.family': 'Lucida Sans Unicode',
          'grid.linestyle': '-',
          'grid.linewidth': 0.5,
          'grid.alpha': 0.75,
          'grid.color': '#ffffff',
          'mathtext.fontset': 'stixsans',
          'mathtext.sf': 'sans',
          'legend.frameon': True,
          'legend.facecolor': '#FFEDCE',
          'figure.dpi': 150}

    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    sns.set_style('darkgrid', rc=rc)
    # colors = sns.color_palette('colorblind', n_colors=8)
    # colors[4] = sns.xkcd_palette(['dusty purple'])[0]
    colors = {'green': '#7AA974', 'light_green': '#BFD598',
              'pale_green': '#DCECCB', 'yellow': '#EAC264',
              'light_yellow': '#F3DAA9', 'pale_yellow': '#FFEDCE',
              'blue': '#738FC1', 'light_blue': '#A9BFE3',
              'pale_blue': '#C9D7EE', 'red': '#D56C55', 'light_red': '#E8B19D',
              'pale_red': '#F1D4C9', 'purple': '#AB85AC',
              'light_purple': '#D4C2D9'}
    if return_colors:
        return colors

def save_seg(fname, image, mask, fill_contours=True, ip_dist=0.160,
             bar_length=10, title=None, colormap='hls'):
    """
    Saves a merge of a segmentation mask and the original image for a
    sanity check.

    Parameters
    ----------
    fname : str
        The file will be saved with this path.
    image : 2d-array, float
        The original image on which the segmentation mask will be overlaid.
    mask : 2d-array, bool
        Boolean segmentation mask of the original image.
    contours: bool
        If True, contours of segmented objects will be filled.
    ip_dist : float
        Interpixel distance for the image. This is used for computing the
        scalebar length.  This should be in units of microns. Default
        value is 0.160 microns per pixel.
    bar_length : int
        The length of the desired scalebar in units of microns.
    title : str, optional
        Title for the image.
    colormap : str
        Colormap for labeling the objects. Default is the high-contrast
        'hls'. This can take any standard colormap string.

    Return
    ------
    fig : Matplotlib Figure object
        Figure containing the axis of the plotted image.
    """

    # Make copies of the image and mask.
    image_copy = np.copy(image)
    mask_copy = np.copy(mask)

    # Burn the scalebar into the upper-left hand  of th image.
    num_pix = int(bar_length / ip_dist)
    image = (image_copy - image_copy.min()) /\
            (image_copy.max() - image_copy.min())
    image[10:20, 10:10 + num_pix] = 1.0

    # Make sure the mask is a boolean image.
    if type(mask) != bool:
        mask = mask_copy > 0

    # Find the contours of the mask.
    conts = skimage.measure.find_contours(mask, 0)

    # Plot the image and generate the contours.
    with sns.axes_style('white'):
        fig = plt.figure()
        plt.imshow(image, cmap=plt.cm.Greys_r)

        # Plot all of the contours
        colors = sns.color_palette(colormap, n_colors=len(conts))
        for i, c in enumerate(conts):
            plt.plot(c[:, 1], c[:, 0], color=colors[i], lw=0.75)
            if fill_contours is True:
                plt.fill(c[:, 1], c[:, 0], color=colors[i], alpha=0.5)

        # Remove the axes.
        plt.xticks([])
        plt.yticks([])

        # Add title if provided.
        if title is not None:
            plt.title(title)

        # Tighten up and save the image.
        plt.tight_layout()
        plt.savefig(fname, bbox_inches='tight')
        plt.close()
    return fig

def show_connections(fname, image, data, title=None, bar_length=10,
                     ip_dist=0.16):
    """
    Saves the original phase contrast image with the segmented
    centroids and the manually recorded markers linked by lines.

    Parameters
    ----------
    fname : str
        Filename to save the image wish shown connections between
        segmented object centroids and the markers.
    image : 2d-array
        Original phase contrast image over which the points will
        be drawn
    data : Pandas DataFrame
        DataFrame containing the marker x and y positions and the
        centroid x and y positions.
    title : str
        Title to be applied to the image. If not specified, none will
        be included.
    bar_length : int
        Length of the scalebar in units of microns. Default value
        is 10.
    ip_dist : float
        Interpixel distance of the image. This should be in units of
        microns per pixel. Default value is 0.16 microns per pixel.

    Returns
    -------
    fig : Matplotlib Figure Canvas
        Figure canvas of the plot.
    """
    # Add the scale bar to the image.
    if image.max() > 1:
        image = (image - image.min()) / (image.max() - image.min())
    num_pix = int(bar_length / ip_dist)
    image_copy = np.copy(image)
    image_copy[10:20, 10:10 + num_pix] = 1.0

    # Define the colors for survivors and corpses.
    colors = {False: '#D56C55', True: '#08AADE'}

    # Group the DataFrame by survival.
    grouped = pd.groupby(data, 'survival')

    # Show the image
    with sns.axes_style('white'):
        fig = plt.figure()
        plt.imshow(image_copy, cmap=plt.cm.Greys_r)
        plt.plot([], [], '-o', ms=3, lw=1, color=colors[True],
                 label='survivor')
        plt.plot([], [], '-o', ms=3, lw=1, color=colors[False], label='goner')
        plt.legend(loc='lower left')
        for g, d in grouped:
            for i in range(len(d)):
                # Parse the positions
                m_x = d.iloc[i]['x_pos']
                m_y = d.iloc[i]['y_pos']
                c_x = d.iloc[i]['label_cent_x']
                c_y = d.iloc[i]['label_cent_y']
                # Plot the connections.
                plt.plot((m_x, c_x), (m_y, c_y), '-', ms=3, lw=1,
                         color=colors[g])
                plt.plot(m_x, m_y, 'o', ms=3, lw=1, color=colors[g])
                plt.plot(c_x, c_y, 'o', ms=3, markerfacecolor='w',
                         markeredgecolor=colors[g], markeredgewidth=1)

        # Format the axes
        plt.xticks([])
        plt.yticks([])

        # Add a title if necessary.
        if title is not None:
            plt.title(title, fontsize=12)
        plt.savefig(fname, bbox_inches='tight')
    return fig

def boilerplate(**kwargs):
    """
    A boiler-plate for a bokeh plotting figure. See `bokeh.plotting.figure`
    for more documentation.
    """
    # Make a bokeh figure axis.
    if kwargs is not None:
        p = bokeh.plotting.figure(**kwargs)
    else:
        p = bokeh.plotting.figure()

    # Apply the styling to the figure axis.
    p.background_fill_color = '#E3DCD0'
    p.grid.grid_line_color = '#FFFFFF'
    p.grid.grid_line_width = 0.75
    p.axis.minor_tick_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.axis_line_color = None
    p.axis.axis_label_text_font = 'Lucida Sans Unicode'
    p.axis.major_label_text_font = 'Lucida Sans Unicode'
    p.axis.axis_label_text_font_style = 'normal'
    p.axis.axis_label_text_font_size = '1em'
    p.axis.major_label_text_font_size = '0.75em'
    p.axis.axis_label_text_color = '#3c3c3c'
    p.axis.axis_label_standoff = 3
    return p



def imshow(im, color_mapper=None, plot_height=400, length_units='pixels',
           interpixel_distance=1.0, return_glyph=False):
    """
    Display an image in a Bokeh figure.
    Parameters
    ----------
    im : 2-dimensional Numpy array
        Intensity image to be displayed.
    color_mapper : bokeh.models.LinearColorMapper instance, default None
        Mapping of intensity to color. Default is 256-level Viridis.
    plot_height : int
        Height of the plot in pixels. The width is scaled so that the
        x and y distance between pixels is the same.
    length_units : str, default 'pixels'
        The units of length in the image.
    interpixel_distance : float, default 1.0
        Interpixel distance in units of `length_units`.
    return_glyph : book, default False
        If True, the image GlyphRenderer is also returned.
    Returns
    -------
    output : bokeh.plotting.figure instance
        Bokeh plot with image displayed.
    Notes
    -----
    We thank Justin Bois for writing this function. http://bebi103.caltech.edu
    """
    # Get shape, dimensions
    n, m = im.shape
    dw = m * interpixel_distance
    dh = n * interpixel_distance

    # Set up figure with appropriate dimensions
    plot_width = int(m / n * plot_height)
    kwargs = {'plot_height': plot_height, 'plot_width': plot_width,
              'x_range': [0, dw], 'y_range': [0, dh],
              'x_axis_label': length_units, 'y_axis_label': length_units}
    p = boilerplate(**kwargs)

    # Set color mapper; we'll do Viridis with 256 levels by default
    if color_mapper is None or color_mapper is 'viridis':
        color_mapper = bokeh.models.LinearColorMapper(
            bokeh.palettes.viridis(256))
    if color_mapper is 'Greys_r':
        color_mapper = bokeh.models.LinearColorMapper(
            bokeh.palettes.gray(256))

    # Display the image
    im_bokeh = p.image(image=[im[::-1, :]], x=0, y=0, dw=dw, dh=dh,
                       color_mapper=color_mapper)

    if return_glyph is True:
        return p, im_bokeh
    else:
        return p


def fill_between(p, domain, range_1, range_2, **kwargs):
    """
    Creates a shaded area between two curves on a bokeh figure
    similar to matplotlibs `fill_between`.
    Parameters
    ----------
    p : Bokeh Figure
        Figure canvas on which to plot shaded region.
    domain : numpy array
        Range of x-values over which to shade
    range_1, range_2 : numpy array
        The two curves to fill between.
    """

    # Ensure that domain, range_1, and range_2 are the same length.
    lengths = np.unique([len(domain), len(range_1), len(range_2)])
    if len(lengths) > 1:
        raise RuntimeError(
            "domain, range_1, and range_2 must be the same length.")

    # Set up the patch objects.
    patch_x = np.append(domain, domain[::-1])
    patch_y = np.append(range_1, range_2[::-1])

    # Add a patch to the bokeh Figure.
    p.patch(patch_x, patch_y, **kwargs)
