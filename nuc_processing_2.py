import numpy as np
import pandas as pd

# skimage stuff
import skimage.feature
import skimage.filters
import skimage.filters.rank
import skimage.io
import skimage.morphology
import skimage.segmentation
from skimage.feature import blob_dog, blob_log, blob_doh
from skimage.color import rgb2gray

import scipy.ndimage

import bebi103

import bokeh

#### HELPER FUNCTIONS ######

def load_img(fname):
    im = skimage.io.imread(fname)[0,:,:]
    return im 

def plot_hist(hist_bin, title, y_axis_type='linear'):
    """Make plot of image histogram."""
    p = bokeh.plotting.figure(plot_height=300,
                              plot_width=400,
                              y_axis_type=y_axis_type,
                              x_axis_label='intensity',
                              y_axis_label='count',
                              title=title)
    hist, bins = hist_bin
    p.line(bins, hist, line_width=2)

    return p

def show_two_ims(im_1, im_2, titles=[None, None], interpixel_distances=[0.13, 0.13],
                 color_mapper=None):
    """Convenient function for showing two images side by side."""
    p_1 = bebi103.viz.imshow(im_1,
                             plot_height=300,
                             title=titles[0],
                             color_mapper=color_mapper,
                             interpixel_distance=interpixel_distances[0],
                             length_units='µm')
    p_2 = bebi103.viz.imshow(im_2,
                             plot_height=300,
                             title=titles[1],
                             color_mapper=color_mapper,
                             interpixel_distance=interpixel_distances[0],
                             length_units='µm')
    p_2.x_range = p_1.x_range
    p_2.y_range = p_1.y_range
    
    return bokeh.layouts.gridplot([p_1, p_2], ncols=2)



#### IMAGE PROCESSING #### 
def granule_filt(im, std = 40.0, selem = skimage.morphology.square(1)):
    im_float = skimage.img_as_float(im)
    im_bg = skimage.filters.gaussian(im_float, std)
    # Subtract background
    im_no_bg = im_float - im_bg
    

    # Use a median filter in a small square structuring element for median filter

    im_snp_filt = skimage.filters.rank.median(im_no_bg, selem)
    return im_no_bg, im_snp_filt


def granule_thresh(im, selem = skimage.morphology.disk(10), white_true=True, k_range=(0.0, 1.0), 
                   min_size=400, selem2 = skimage.morphology.disk(1)):
    """
    Threshold image as described above.  Morphological mean filter is 
    applied using selem.
    """    
    # Determine comparison operator
    if white_true:
        compare = np.greater
        sign = -1
    else:
        compare = np.less
        sign = 1
    
    # Do the mean filter
    im_mean = skimage.filters.rank.mean(im, selem)

    # Compute number of pixels in binary image as a function of k
    k = np.linspace(k_range[0], k_range[1], 100)
    n_pix = np.empty_like(k)
    for i in range(len(k)):
        n_pix[i] = compare(im, k[i] * im_mean).sum() 

    # Compute rough second derivative
    dn_pix_dk2 = np.diff(np.diff(n_pix))

    # Find index of maximal second derivative
    max_ind = np.argmax( dn_pix_dk2)
    
    # Use this index to set k
    k_opt = k[max_ind - sign * 2]
    print('k_opt: ', k_opt)
 
    # Threshold with this k
    im_bw = compare(im, k_opt * im_mean)
    # Remove all the small objects
    im_bw = skimage.morphology.remove_small_objects(im_bw, min_size=min_size)
    
    # close binary image
    im_bw = skimage.morphology.binary_closing(im_bw, skimage.morphology.disk(10))



    return im_bw, k_opt

def get_blobs(im, min_sigma = 1, max_sigma = 2, threshold = 0.001):
    image_gray = rgb2gray(im)
    blobs_log = blob_dog(image_gray, min_sigma = min_sigma , max_sigma=max_sigma , 
                         threshold=threshold)
    blobs_log[:, 2] = blobs_log[:, 2] * np.sqrt(2)
    p = np.zeros(im.shape)
    for coords in blobs_log:
 
        r, c, rad = coords.astype(int)
  
        p[r-rad:r+rad, c-rad:c+rad] = 1
    return blobs_log, p 


def get_nuclei(im_bw, im, im_black, blobs = False, threshR = 1., threshS = 1800, std = 68.0, min_sigma = 1,
               max_sigma = 2, threshold = 0.001, background = True ):
    # Label binary image; backward kwarg says value in im_bw to consider backgr.
    im_labeled, n_labels = skimage.measure.label(
                            im_bw, background=0, return_num=True)
    im_float = skimage.img_as_float(im)
    im_bg = skimage.filters.gaussian(im_float, 100)
    # Subtract background
    
    im_no_bg = im_float 

    im_no_bg = im_float - im_bg

    print('Number of individual objects = ', n_labels - 1)

    
    # Get properties about the YFP channel
    im_props = skimage.measure.regionprops(im_labeled, intensity_image=im_black)
    data = [[prop.label, prop.mean_intensity] for prop in im_props]

    df = pd.DataFrame(data=data, columns=['label', 'mean intensity (a.u.)'])
    im_filt = np.copy(im_labeled)
    bad = np.where(df['mean intensity (a.u.)'] > 1200)[0]

    for i in bad:
        
        inds = np.where(im_labeled == i + 1, im_labeled, 0)
        im_filt = im_filt -  inds
    # background
    im_props = skimage.measure.regionprops(im_labeled, intensity_image=im_no_bg)
    data = [[prop.label, prop.area, prop.mean_intensity, prop.moments_hu[1], 
             prop.bbox[0], prop.bbox[1]] for prop in im_props]
    df = pd.DataFrame(data=data, columns=['label', 'area (sq ip distance)', 
                                          'mean intensity (a.u.)', 'nuclei', 'bbox_min', 'bbox_max'])
    df = df.loc[np.unique(im_filt)[1:]-1, :]

    if blobs:
        
        blobs_log, p = get_blobs(im, min_sigma = min_sigma , max_sigma=max_sigma , 
                         threshold=threshold)
        
        im_bw_blobs = im_bw*p
        
        im_labeled_blobs, _ = skimage.measure.label(
                            im_bw_blobs, background=0, return_num=True)
        im_props_blobs = skimage.measure.regionprops(im_labeled_blobs, intensity_image=im_no_bg)
        
        data = [[prop.label, prop.area, prop.mean_intensity, prop.moments_hu[1], 
             prop.bbox[0], prop.bbox[1]] for prop in im_props_blobs]
        
        df_blobs = pd.DataFrame(data=data, columns=['label', 'area (sq ip distance)', 
                                          'mean intensity (a.u.)', 'nuclei', 'bbox_min', 'bbox_max']) 
        df_blobs['nuclei id'] = 0
        nuc_blobs = np.zeros(len(df_blobs))
        for nuc in df['label'].values:
            inds = np.where(im_filt == nuc, im_filt, 0)/nuc
            blobs_inds = np.unique(im_labeled_blobs*inds)[1:]
            
            nuc_blobs[blobs_inds.astype(int) - 1] = nuc
        df_blobs['nuclei id'] = nuc_blobs
        im_filt_blobs = im_filt*im_bw_blobs
        df_blobs = df_blobs.loc[df_blobs['nuclei id'] > 0, :]
        return im_filt, im_filt_blobs, df_blobs
    return im_labeled, im_filt, df


def full_process(fname, selem = skimage.morphology.disk(10), selem2 = skimage.morphology.disk(2), 
                 std = 40.0, min_size = 400, threshR = 1., threshS = 1800, blobs = False, background = True):
    im = load_img(fname)

    try:
        im_black = load_img(fname[:-4] + '_black.tif')
    except: 
        im_black = skimage.io.imread(fname[:-4] + '_black.tif')
    _, im_snp_filt = granule_filt(im, std = std)
    im_bw, _ = granule_thresh(im_snp_filt, selem=selem,
                   min_size=min_size, selem2 = selem2)
    im_labeled, im_filt, df = get_nuclei(im_bw, im, im_black, threshR = threshR, threshS = threshS, blobs = blobs,
                                        background = True )
    return im_labeled, im_filt, im, df
    

    
    
    