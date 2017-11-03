# -*- coding: utf-8 -*-
import numpy as np
import skimage.io
import skimage.segmentation
import skimage.morphology
import skimage.measure
import skimage.filters
import scipy.ndimage
import xmltodict
import json
import glob
import pandas as pd


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

# #################


def find_zero_crossings(im, selem, thresh):
    """
    This  function computes the gradients in pixel values of an image after
    applying a sobel filter to a given image. This  function is later used in
    the Laplacian of Gaussian cell segmenter (log_segmentation) function. The
    arguments are as follows.

    Parameters
    ----------
    im : 2d-array
        Image to be filtered.
    selem : 2d-array, bool
        Structural element used to compute gradients.
    thresh :  float
        Threshold to define gradients.

    Returns
    -------
    zero_cross : 2d-array
        Image with identified zero-crossings.

    Notes
    -----
    This function as well as `log_segmentation` were written by Justin Bois.
    http://bebi103.caltech.edu/
    """

    # apply a maximum and minimum filter to the image.
    im_max = scipy.ndimage.filters.maximum_filter(im, footprint=selem)
    im_min = scipy.ndimage.filters.minimum_filter(im, footprint=selem)

    # Compute the gradients using a sobel filter.
    im_filt = skimage.filters.sobel(im)

    # Find the zero crossings.
    zero_cross = (((im >= 0) & (im_min < 0)) | ((im <= 0) & (im_max > 0)))\
        & (im_filt >= thresh)

    return zero_cross


# #################
def log_segmentation(im, selem=None, thresh=0.0001, radius=2.0,
                     median_filt=True, clear_border=True, label=False):
    """
    This function computes the Laplacian of a gaussian filtered image and
    detects object edges as regions which cross zero in the derivative.

    Parameters
    ----------
    im :  2d-array
        Image to be processed. Must be a single channel image.
    selem : 2d-array, bool
        Structural element for identifying zero crossings. Default value is
        a 2x2 pixel square.
    radius : float
        Radius for gaussian filter prior to computation of derivatives.
    median_filt : bool
        If True, the input image will be median filtered with a 3x3 structural
        element prior to segmentation.
    selem : 2d-array, bool
        Structural element to be applied for laplacian calculation.
    thresh : float
        Threshold past which
    clear_border : bool
        If True, segmented objects touching the border will be removed.
        Default is True.
    label : bool
        If True, segmented objecs will be labeled. Default is False.

    Returns
    -------
    im_final : 2d-array
        Final segmentation mask. If label==True, the output will be a integer
        labeled image. If label==False, the output will be a bool.

    Notes
    -----
    We thank Justin Bois in his help writing this function.
    https://bebi103.caltech.edu
    """

    # Test that the provided image is only 2-d.
    if len(np.shape(im)) > 2:
        raise ValueError('image must be a single channel!')

    # Determine if the image should be median filtered.
    if median_filt is True:
        selem = skimage.morphology.square(3)
        im_filt = scipy.ndimage.median_filter(im, footprint=selem)
    else:
        im_filt = im
    # Ensure that the provided image is a float.
    if np.max(im) > 1.0:
        im_float = skimage.img_as_float(im_filt)
    else:
        im_float = im_filt

    # Compute the LoG filter of the image.
    im_LoG = scipy.ndimage.filters.gaussian_laplace(im_float, radius)

    # Define the structural element.
    if selem is None:
        selem = skimage.morphology.square(3)

    # Using find_zero_crossings, identify the edges of objects.
    edges = find_zero_crossings(im_LoG, selem, thresh)

    # Skeletonize the edges to a line with a single pixel width.
    skel_im = skimage.morphology.skeletonize(edges)

    # Fill the holes to generate binary image.
    im_fill = scipy.ndimage.morphology.binary_fill_holes(skel_im)

    # Remove small objects and objects touching border.
    im_final = skimage.morphology.remove_small_objects(im_fill)
    if clear_border is True:
        im_final = skimage.segmentation.clear_border(im_final, buffer_size=5)

    # Determine if the objects should be labeled.
    if label is True:
        im_final = skimage.measure.label(im_final)

    # Return the labeled image.
    return im_final


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
