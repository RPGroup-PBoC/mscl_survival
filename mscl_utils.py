import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import mscl_utils as mscl

import skimage.io
import skimage.measure
import skimage.morphology
import skimage.segmentation
import skimage.filters
import scipy.ndimage
import os

import bokeh.plotting
import cairosvg

import pandas as pd
import xmltodict
import json


# ----------------------------------------------------------------------------
# Image Processing and Marker Linking Utilities
# ----------------------------------------------------------------------------

def compute_mean_bg(phase_image, fluo_image, method='isodata', obj_dark=True):
    """
    Computes the mean background fluorescence of the inverted segmentation
    mask.

    Parameters
    ----------
    phase_image : 2d-array, int or float.
        The phase contrast image used for generating the inverse segmentation
        mask. If this image is not a float with pixel values in (0, 1), it
        will be renormalized.
    fluo_image : 2d-array, int
        The fluorescence image used to calculate the mean pixel value. If
        flatfield correction is necessary, it should be done before this
        sending to this function.
    method: string, ['otsu', 'yen', 'li', 'isodata'], default 'isodata'
        Automated thresholding method to use. Default is 'isodata' method.
    obj_dark : bool, default True
        If True, objects will be **darker** than the automatically generated
        threshold value. If False, objects are deemed to be brighter.

    Returns
    -------
    mean_bg: float
        The mean background fluorescence of the image.
    """

    # Ensure that the image is renormalized.
    if (phase_image > 1.0).any():
        phase_image = (phase_image - phase_image.min()) /\
                      (phase_image.max() - phase_image.min())
    # Perform the background subtraction.
    im_blur = skimage.filters.gaussian(phase_image, sigma=50)
    im_sub = phase_image - im_blur

    # Determine the method to use.
    methods = {'otsu': skimage.filters.threshold_otsu,
               'yen': skimage.filters.threshold_yen,
               'li': skimage.filters.threshold_li,
               'isodata': skimage.filters.threshold_isodata}

    # Determine the threshold value.
    thresh_val = methods[method](im_sub)

    # Generate the inverted segmentation mask and dilate.
    if obj_dark is True:
        im_thresh = im_sub < thresh_val
    else:
        im_thresh = im_sub > thresh_val

    selem = skimage.morphology.disk(20)
    im_dil = skimage.morphology.dilation(im_thresh, selem=selem)

    # Mask onto the fluroescence image and compute the mean background value.
    mean_bg = np.mean(fluo_image[im_dil < 1])
    return mean_bg


def median_flatfield(image_stack, medfilter=True, selem='default',
                     return_profile=False):
    """
    Computes a illumination profile from the median of all images
    and corrects each individual image.

    Parameters
    ----------
    image_stack: scikit-image ImageCollection
        Series of images to correct. The illumination profile is created
        from computing the median filter of all images in this collection.
    medfilter: bool, default True
        If True, each individiual image will be prefiltered using a median
        filter with  a given selem.
    selem : string or structure, default 3x3 square
        Structural element to use for the median filtering. Default  is
        a 3x3 pixel square.
    return_profile: bool, default False
        If True, the illumination profiled image will be returned.

    Returns
    -------
    ff_ims : list of 2d-array
        Flatfield corrected images.
    med_im : 2d-array
        Illumination profile produced from the median of all images in
        image stack.
    """

    # Determine if the prefiltering should be performed.
    if medfilter is True:

        # Define the structural element.
        if selem is 'default':
            selem = skimage.morphology.square(3)
        image_stack = [scipy.ndimage.median_filter(
            im, footprint=selem) for im in image_stack]

    # Compute the median filtered image.
    med_im = np.median(image_stack, axis=0)

    # Perform the correction.
    ff_ims = [(i / med_im) * np.mean(med_im) for i in image_stack]

    if return_profile is True:
        return [ff_ims, med_im]
    else:
        return ff_ims


def average_stack(im, median_filt=True):
    """
    Computes an average image from a provided array of images.

    Parameters
    ----------
    im : list or arrays of 2d-arrays
        Stack of images to be filtered.
    median_filt : bool
        If True, each image will be median filtered before averaging.
        Median filtering is performed using a 3x3 square structural element.

    Returns
    -------
    im_avg : 2d-array
        averaged image with a type of int.
    """

    # Determine if the images should be median filtered.
    if median_filt is True:
        selem = skimage.morphology.square(3)
        im_filt = [scipy.ndimage.median_filter(i, footprint=selem) for i in im]
    else:
        im = im_filt

    # Generate and empty image to store the averaged image.
    im_avg = np.zeros_like(im[0]).astype(int)
    for i in im:
        im_avg += i
    im_avg = im_avg / len(im)
    return im_avg


def generate_flatfield(im, im_field, median_filt=True):
    """
    Corrects illumination of a given image using a dark image and an image of
    the flat illumination.

    Parameters
    ----------
    im : 2d-array
        Image to be flattened.
    im_field: 2d-array
        Average image of fluorescence illumination.
    median_filt : bool
        If True, the image to be corrected will be median filtered with a
        3x3 square structural element.

    Returns
    -------
    im_flat : 2d-array
        Image corrected for uneven fluorescence illumination. This is performed
        as

        im_flat = (im  / im_field ) * mean(im_field)

    Raises
    ------
    RuntimeError
        Thrown if bright image and dark image are approximately equal. This
        will result in a division by zero.
    """
    # Compute the mean field value.
    mean_diff = np.mean(im_field)

    if median_filt is True:
        selem = skimage.morphology.square(3)
        im_filt = scipy.ndimage.median_filter(im, footprint=selem)
    else:
        im_filt = im

    # Compute and return the flattened image.
    im_flat = (im_filt / im_field) * mean_diff
    return im_flat

# for segmentation.


def contour_seg(image, level=0.3, selem='default', perim_bounds=(5, 1E3),
                ip_dist=0.160, ecc_bounds=(0.7, 1), area_bounds=(1, 50),
                return_conts=False, min_int=0.2):
    """
    Identifies contours around dark objects in a phase contrast image.

    Parameters
    ----------
    image: 2d-array
        Phase contrast image of interest.
    level: float
        Level at which to draw contours on black top-hat filtered image.
        Default value is 0.3.
    selem: 2d-array or string
        Structuring element to use for the black top-hat filtering procedure
        Default value is a disk with a diameter of 20 pixels.
    perim_bounds: length 2 tuple
        Lower and upper perimeter bounds of approved objects. This should be
        in units of microns. The default values are 5 and 25 microns for the
        lower and upper bound, respectively.
    ip_dist : float
        Interpixel distance of the image in units of microns per pixel. The
        default value is 0.160 microns per pixel.
    area_bounds : tuple of float
        Upper and lower bounds for selected object areas. These should be
        given in units of square microns.
    ecc_bounds : tuple of float
        Bounds for object eccentricity. Default values are between 0.5 and 1.0.
    return_conts : bool
        If True, the x and y coordinates of the individual contours will be
        returned. Default value is False

    Returns
    -------
    im_lab : 2d-array, int
        Two dimensional image where each individual object is labeled.

    conts : 1d-array
        List of contour coordinates. Each entry of this array comes as
        an x,y pair of arrays. Has the same length as the number of
        contoured objects. This is only returned if `return_conts` is
        True.

    """

    # Apply the white top-hat filter.
    if selem == 'default':
        selem = skimage.morphology.disk(20)

    # Normalize the image.
    image = (image - image.min()) / (image.max() - image.min())

    # Blur and background subtract the image.
    im_blur = skimage.filters.gaussian(image, sigma=5)
    im_sub = image - im_blur

    # Apply the black tophat filter.
    im_filt = skimage.morphology.black_tophat(im_sub, selem)

    # Find the contours and return.
    conts = skimage.measure.find_contours(im_filt, level)

    # Make an empty image for adding the approved objects.
    objs = np.zeros_like(image)

    # Loop through each contour.
    for _, c in enumerate(conts):
        perim = 0
        for j in range(len(c) - 1):
            # Compute the distance between points.
            distance = np.sqrt((c[j + 1, 0] - c[j, 0])**2 +
                               (c[j + 1, 1] - c[j, 1])**2)
            perim += distance * ip_dist

        # Test if the perimeter is allowed by the user defined bounds.
        if (perim > perim_bounds[0]) & (perim < perim_bounds[1]):

            # Round the contours.
            c_int = np.round(c).astype(int)

            # Color the image with the contours and fill.
            objs[c_int[:, 0], c_int[:, 1]] = 1.0

    # Fill and label the objects.
    objs_fill = scipy.ndimage.binary_fill_holes(objs)
    objs_fill = skimage.morphology.remove_small_objects(objs_fill)
    im_lab = skimage.measure.label(objs_fill)

    # Apply filters.
    approved_obj = np.zeros_like(im_lab)
    props = skimage.measure.regionprops(im_lab, image)
    for prop in props:
        area = prop.area * ip_dist**2
        ecc = prop.eccentricity
        if (area < area_bounds[1]) & (area > area_bounds[0]) &\
            (ecc < ecc_bounds[1]) & (ecc > ecc_bounds[0]) &\
                (prop.mean_intensity < min_int):
            approved_obj += (im_lab == prop.label)
    im_lab = skimage.measure.label(approved_obj)

    if return_conts is True:
        return conts, im_lab
    else:
        return im_lab


def marker_parse(fname, type_dict={1: False, 2: True}):
    """
    Parses the XML file produced from the CellCounter ImageJ plugin and
    packages the marker positions and type into a Pandas DataFrame.

    Parameters
    ----------
    fname : str
        Path to the XML file of interest.
    type_dict : dict
        Dictionary of types and survival. Default is assigning type 1
        as death and type 2 as survival.

    Returns
    -------
    df : Pandas DataFrame
        Data frame containing x and y positions of markers as well as
        the type classification.
    """
    with open(fname, 'r') as f:
        positions = xmltodict.parse(f.read())

    # Extract only the marker data.
    markers = positions['CellCounter_Marker_File']['Marker_Data']['Marker_Type']

    # Find the total number of types and loop through them to make data frames.
    dfs = []
    num_types = len(markers)
    for i in range(num_types):
        try:
            type_marks = markers[i]['Marker']
            _df = pd.DataFrame(type_marks)
            # Insert a column keeping track of the type
            _df.insert(0, 'survival', type_dict[int(markers[i]['Type'])])
            dfs.append(_df)
        except:
            pass

    # Concatenate the data frames ignorning indexing.
    df = pd.concat(dfs, axis=0, ignore_index=True)

    # Clean up the data frame and return.
    df.drop('MarkerZ', 1, inplace=True)
    df.columns = ['survival', 'x_pos', 'y_pos']
    df['x_pos'] = df['x_pos'].astype(int)
    df['y_pos'] = df['y_pos'].astype(int)
    return df


def link_markers(markers, seg_mask, fluo_image, ip_dist=0.160,
                 return_coords=False, max_dist=5,
                 position_labels=('x_pos', 'y_pos')):
    """
    Maps markers from one image to centroids of segmented objects from
    another. This assumes a marker belongs to the object with the minimum
    marker-centroid distance.

    Parameters
    ----------
    markers : Pandas DataFrame
        DataFrame containing the x and y positions of the markers.
    seg_mask : 2d-array, int
        Labeled segmentation mask. The coordinates of the object centroids
        will be calculated from this image.
    fluo_image : 2d-array, float or int
        The fluorescence image used to extract intensities. If None, no
        intensity information will be returned. These intensity values
        will be returned as an intensity per square physical distance
        as given by `ip_dist`.
    ip_dist :  float
        Interpixel distance for the image. Default value is 0.160 microns
        per pixel.
    return_coords : bool
        If True, the paired coordinates will be returned as a tuple. It
        will have the form ((mark_x, mark_y), (cent_x, cent_y)). Default
        value is False.
    max_dist : float
        Maximum distance to keep. Default Value is 5 microns.
    position_labels :  tuple of str
        Labels of position markers in the markers DataFrame in the order
        of x position and y position. Default is `x_pos` and `y_pos`.

   Returns
   -------
   df : Pandas DataFrame
       DataFrame containing survival type, marker positions, mask label,
       area, and intensity if provided. Note this is not returned if
       `in_place = True`.
   coords : list of tuple
       A list of tuples containing the marker x,y positions and the
       coordinates of the associated segmentation centroid. This
       is only returned if `return_coords`==True.

    """
    # Compute the properties from the segmentation mask.
    props = skimage.measure.regionprops(seg_mask, fluo_image)
    area, intensity, labels, centroids = [], [], [], []
    for prop in props:
        area.append(prop.area * ip_dist**2)
        intensity.append(prop.mean_intensity / ip_dist**2)
        labels.append(prop.label)
        centroids.append(prop.centroid)

    # Set up a list to store the coordinates and duplicate the df.
    coords = []
    if type(markers) == str:
        df = pd.DataFrame([intensity, area]).T
        df.columns = ['intensity', 'area']
        df.insert(np.shape(df)[1], 'dist', 0)
        df.insert(0, 'label_cent_y', 0)
        df.insert(0, 'label_cent_x', 0)
        df.insert(0, 'mask_label', labels)
        df.insert(0, 'y_pos', 0)
        df.insert(0, 'x_pos', 0)
        df.insert(0, 'survival', False)
        return df

    else:
        df = markers.copy(deep=True)

        # Compute the minimum distances.
        for i in range(len(markers)):
            distances = []
            x = markers.iloc[i][position_labels[0]]
            y = markers.iloc[i][position_labels[1]]

            # Loop through each centroid and find the minimum distance.
            for c in centroids:
                dist = np.sqrt((x - c[1])**2 + (y - c[0])**2)
                distances.append(dist)
            if len(distances) == 0:
                df.set_value(i, 'dist', 1E6)
                pass
            else:
                # Find the index with the minimum distance.
                min_ind = np.argmin(distances)
                coords.append(((x, y),
                               (centroids[min_ind][1], centroids[min_ind][0])))

                # Determine if a new DataFrame should be made or not.
                # TODO: Make do this by generating a dictionary and appending
                # to an existing dataframe.
                df.set_value(i, 'mask_label', labels[min_ind])
                df.set_value(i, 'label_cent_x', centroids[min_ind][1])
                df.set_value(i, 'label_cent_y', centroids[min_ind][0])
                df.set_value(i, 'intensity', intensity[min_ind])
                df.set_value(i, 'area', area[min_ind])
                df.set_value(i, 'dist', distances[min_ind])

        # Apply the distance filter.
        df = df[df['dist'] <= (max_dist / ip_dist)]
        if return_coords is True:
            return df, coords
        else:
            return df


def scrape_metadata(fname, channels=('Brightfield', 'GFP'), return_date=True):
    """
    Takes an image metadata file and returns the datea nd GFP exposure time.

    Parameters
    ----------
    fname : str
        Pat of the metadata file to parse.
    channel : tuple of str
        The channels from which to scrape the exposure time. A single channel
        name can be given. Default is ('Brightfield', 'GFP').
    return_date : bool
        If True, the date of the acquisition will also be returned.

    Returns
    -------
    exposure: dict or float
        The exposure time of the desired channel. If multiple channels are
        given, this will be a tuple of the exposure times. If return_date is
        True, the date will also be in this dictionary.
   """

    # Open the metadata file.
    with open(fname, 'r') as f:
        metadata = json.load(f)

    # Get a list of the keys in the metadata file.
    keys = metadata.keys()

    # Determine if a single channel or multiple channel exposures are desired.
    if (type(channels) != tuple) & (type(channels) != str):
        raise TypeError('desired channels must be a tuple or a string.')
    else:
        if type(channels) == str:
            num_channels = 1
            channels = (channels)
            exposure = None
        else:
            num_channels = len(channels)
            exposure = []

    # Loop through each desired channel and scrape the exposure.
    for i in range(num_channels):
        for k in keys:
            try:
                chan = metadata[k]['Channel']
                if chan.lower() == channels[i].lower():
                    _exposure = metadata[k]['Exposure-ms']
                if num_channels == 1:
                    exposure = _exposure
                else:
                    if i == 0:
                        exposure = {channels[i] + '_exp_ms': _exposure}
                    else:
                        exposure[channels[i] + '_exp_ms'] = _exposure
            except:
                pass

    if return_date is True:
        # Get the date from the Summary field.
        date = metadata['Summary']['Date'].split('-')
        date = ''.join(date)
        exposure['date'] = date
    return exposure

# ---------------------------------------------------------------------------
# Plotting utilities using Bokeh and Matplotlib
# ---------------------------------------------------------------------------


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

# TODO: Rewrite this with a  Bokeh backend.


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


# Plotting utilities
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
          #   'xtick.labelsize': 'large',
          #   'ytick.labelsize': 'large',
          'font.family': 'Lucida Sans Unicode',
          'grid.linestyle': '-',
          'grid.linewidth': 0.5,
          'grid.alpha': 0.75,
          'grid.color': '#ffffff',
          'mathtext.fontset': 'stixsans',
          'mathtext.sf': 'sans',
          'legend.frameon': True,
          'legend.facecolor': '#FFEDCE',
          #   'figure.figsize': (8, 6),
          'figure.dpi': 300}

    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    sns.set_style('darkgrid', rc=rc)
    sns.set_context('paper', rc=rc)
    colors = sns.color_palette('colorblind', n_colors=8)
    colors[4] = sns.xkcd_palette(['dusty purple'])[0]

    if return_colors:
        return colors


# For Bokeh Styling.
def bokeh_boiler(**kwargs):
    # Make a bokeh figure axis.
    if kwargs is not None:
        p = bokeh.plotting.figure(**kwargs)
    else:
        p = bokeh.plotting.figure()

    # Apply the styling to the figure axis.
    p.background_fill_color = '#E3DCD0'
    p.grid.grid_line_color = '#FFFFFF'
    # p.grid.grid_line_dash = 'dotted'
    p.grid.grid_line_width = 0.75
    p.axis.minor_tick_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.axis_line_color = None
    p.axis.axis_label_text_font = 'Lucida Sans Unicode'
    p.axis.major_label_text_font = 'Lucida Sans Unicode'
    p.axis.axis_label_text_font_style = 'normal'
    p.axis.axis_label_text_font_size = '13pt'
    p.axis.major_label_text_font_size = '10pt'
    p.axis.axis_label_text_color = '#3c3c3c'
    p.axis.axis_label_standoff = 3
    p.output_backend = 'svg'
    return p


def bokeh_to_pdf(p, fname):
    bokeh.io.export_svgs(p, '._tmp.svg')
    cairosvg.svg2pdf(url='./._tmp.svg', write_to=fname)
    os.remove('./._tmp.svg')
    print('Saved bokeh figure as {0}'.format(fname))


# For bokeh image display. Shamelessly taken from Justin Bois.
def bokeh_imshow(im, color_mapper=None, plot_height=400, length_units='pixels',
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
    """
    # Get shape, dimensions
    n, m = im.shape
    dw = m * interpixel_distance
    dh = n * interpixel_distance

    # Set up figure with appropriate dimensions
    plot_width = int(m / n * plot_height)
    kwargs = {'plot_height': plot_height, 'plot_width': plot_width,
              'x_range': [0, dw], 'y_range': [0, dh],
              'x_axis_label': length_units, 'y_axis_label': length_units,
              'tools': 'pan, box_zoom, wheel_zoom, reset, resize'}
    p = bokeh_boiler(**kwargs)

    # Set color mapper; we'll do Viridis with 256 levels by default
    if color_mapper is None:
        color_mapper = bokeh.models.LinearColorMapper(
            bokeh.palettes.viridis(256))

    # Display the image
    im_bokeh = p.image(image=[im[::-1, :]], x=0, y=0, dw=dw, dh=dh,
                       color_mapper=color_mapper)

    if return_glyph is True:
        return p, im_bokeh
    else:
        return p


# ---------------------------------------------------------------------------
# MCMC and Other Inferencial Utilities
# ---------------------------------------------------------------------------
def ecdf(data):
    """
    Computes the empirical cumulative distribution function of a data set.

    Parameters
    ----------
    data: pandas Series, Slice, or 1d-array.
        Data set from which the ECDF will be computed. This must be
        one-dimensional

    Returns
    -------
    x, y : 1d-arrays
        x is the sorted values.
        y is the fractional representation of each value.
    """

    x, y = np.sort(data), np.arange(0, len(data), 1) / len(data)
    return (x, y)


def density_binning(data, groupby='shock_group', channel_bin=10, min_cells=20,
                    channel_key='channel_density', survival_key='survival'):
    """
    Bins survival data by a given channel density.

    Parameters
    ----------
    data : pandas DataFrame
        DataFrame containing data with computed channel density,
        survival classifier, and shock speed designation.
    groupby : list of strings.
        Keys by which to group the survival data. Default is 'shock_group'
    channel_bin : float or int
        Bin width for channel density. Default is 10 channels per unit area.
    min_cells : int
        Minimum number of cells to consider for each bin.
    channel_key : string
        Column name for channel density. Default is 'channel_density'.
    survival_key : string
        Column name for survival identifier. Default is 'survival'.

    Returns:
    --------
    bin_data : pandas DataFrame
        Data frame with binned data.
    """

    # Set the bounds for the bins.
    lower_bound = 0
    upper_bound = int(data[channel_key].max())
    bins = np.arange(lower_bound, upper_bound + channel_bin, channel_bin)

    # Sort the data by channel density
    sorted_data = data.sort_values(by=channel_key)

    # Partition into the bins.
    sorted_data['bin_number'] = 0
    sorted_data.reset_index(inplace=True)
    for i in range(1, len(bins) - 1):
        # Assign bin numbers based on channel density
        inds = (sorted_data['channel_density'] >= bins[i - 1]
                ) & (sorted_data['channel_density'] < bins[i + 1])
        sorted_data.loc[inds,  'bin_number'] = i

    # Ensure that the bin numbering scheme is sequential.
    bin_data = sorted_data.copy()
    grouped = bin_data.groupby(groupby)
    for g, d in grouped:
        seq_change = {}
        bin_nos = d['bin_number'].unique()
        for i, b in enumerate(bin_nos):
            bin_data.loc[(bin_data['bin_number'] == b) &
                         (bin_data[groupby] == g), 'bin_number'] = i

    # Regroup the data and congeal bins to a minimum cell number.
    grouped = bin_data.groupby(groupby)
    for g, d in grouped:
        # Group by bin number.
        _grouped = d.groupby('bin_number')[survival_key].count()

        # Find those less than the minimum cell number.
        bin_counts = _grouped.to_dict()
        low_bins = _grouped[_grouped < min_cells].to_dict()

        # Get just the bin numbers.
        bin_nos = list(low_bins.keys())

        # Identify the edges of sequential bins with low cell counts.
        sequential = np.where(np.diff(bin_nos) > 1)[0]
        if (len(sequential) == 0) & (len(bin_nos) != 0):
            paired = [bin_nos]
        else:
            # Need to do fancy indexing here so it returns even single bins.
            paired = [bin_nos[:sequential[0] + 1]]
            _paired = ([bin_nos[sequential[j - 1] + 1:sequential[j] + 1]
                        for j in range(1, len(sequential))])
            for _p in _paired:
                paired.append(_p)
            paired.append(bin_nos[sequential[-1] + 1:])

        # Loop through each pair and determine if they can meet the minimum.
        change_bins = {}

        for i, pair in enumerate(paired):
            if len(pair) > 1:
                summed = np.sum([bin_counts[p] for p in pair])
                if summed >= min_cells:
                    for z in pair:
                        change_bins[z] = pair[0]
                else:
                    for z in pair:
                        change_bins[z] = pair[0] - 1
            else:
                # Deal with edge cases of first and last bin.
                if pair[0] == 1:
                    change_bins[pair[0]] = pair[0] + 1
                elif pair[0] == sorted_data['bin_number'].max():
                    change_bins[pair[0]] = pair[0] - 1

        # Loop through the changed bins and change the value of the bin number
        # in the original dataframe.
        keys = change_bins.keys()
        for key in keys:
            bin_data.loc[(bin_data[groupby] == g) &
                         (bin_data['bin_number'] == key),
                         'bin_number'] = change_bins[key]
    return bin_data


def compute_survival_stats(df, groupby=['shock_group', 'bin_number']):
    """
    Computes the statistics of survival probabilitiy, number of cells, and
    binomial error given a dataframe with binned events. This should be used
    as an apply function on a pandas groupby method.
    """
    def binomial_probability(df):
        n = np.sum(df == True)
        N = len(df)
        return n / N

    def binomial_err(df):
        n = np.sum(df == True)
        N = len(df)
        return np.sqrt(n * (N - n) / N**3)

    def _compute_stats(df):
        stats_dict = dict(prob=binomial_probability(df['survival']),
                          err=binomial_err(df['survival']),
                          mean_chan=df['channel_density'].mean(),
                          n_cells=len(df), n_suv=np.sum(df['survival']))
        return pd.Series(stats_dict)
    grouped = df.groupby(groupby).apply(_compute_stats)
    return pd.DataFrame(grouped).reset_index()


def trace_to_df(trace, model):
    """
    Converts the trace from a pymc3 sampler to a
    Pandas DataFrame.
    """
    def compute_logp(chain):
        """
        Computes the log probability of the provided trace
        at a given chain.
        """
        names = trace.varnames
        var_dict = {}
        for n in names:
            var_dict[n] = trace.get_values(n, chains=chain)
        sample_df = pd.DataFrame(var_dict)

        logp = [model.logp(sample_df.iloc[step]
                           ) for step in range(len(sample_df))]
        return logp

    chains = trace.chains
    for c in tqdm(chains, desc='Processing chains'):
        logp = compute_logp(c)
        if c == 0:
            df = pm.trace_to_dataframe(trace, chains=c)
            df.insert(np.shape(df)[1], 'logp', logp)
        else:
            _df = pm.trace_to_dataframe(trace, chains=c)
            _df.insert(np.shape(_df)[1], 'logp', logp)
            df.append(_df, ignore_index=True)

    return df


def compute_mcmc_statistics(df, ignore_vars='logp'):
    """
    Computes the mode and highest probability density (hpd)
    of the parameters in a given dataframe.
    """
    # Set up the multi indexing.
    var_names = np.array(df.keys())
    if ignore_vars is not None:
        var_names = var_names[var_names != ignore_vars]

    # Generate arrays for indexing and zip as tuples.
    names = [var for var in var_names] * 3
    stats = ['mode', 'hpd_min', 'hpd_max']
    stats = np.array([[s] * len(var_names) for s in stats]).flatten()
    tuples = list(zip(*[names, stats]))

    # Define the index.
    index = pd.MultiIndex.from_tuples(tuples, names=['var', 'stat'])

    # Determine the mode for each
    mode_ind = np.argmax(df['logp'])
    stat_vals = [df.iloc[mode_ind][var] for var in var_names]
    # Compute the min and max vals of the HPD.
    hpd_min, hpd_max = [], []
    for i, var in enumerate(var_names):
        _min, _max = hpd(df[var], 0.95)
        hpd_min.append(_min)
        hpd_max.append(_max)
    for _ in hpd_min:
        stat_vals.append(_)
    for _ in hpd_max:
        stat_vals.append(_)

    # Add them to the array for the multiindex
    flat_vals = np.array([stat_vals]).flatten()
    var_stats = pd.Series(flat_vals, index=index)
    return var_stats
