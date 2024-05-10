from skimage.metrics import structural_similarity as compare_ssim
from skimage.metrics import structural_similarity as ssim
from datetime import  date, datetime,timedelta
from astropy.coordinates import SkyCoord
from skimage.transform import rescale
from sunpy.net import Fido, attrs as SAF
from  scipy.stats.stats import pearsonr

from astropy.time import Time
from sunpy.net import Fido
from sunpy.net import attrs as a
import sunpy.map

import matplotlib
from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np 



import scipy.ndimage as scnd
import sunpy.map
import numpy as np 
import sunpy
import glob 
import os
import warnings
warnings.filterwarnings("ignore")



#-----------------------------------------------------------------------------------------
#------------------------------------- AUX Functions -------------------------------------
#-----------------------------------------------------------------------------------------
def creating_time_vector(iris_file):
    """
    Create the time vector in datetime.strptime format for IRIS data cubes.

    Input:
    - iris_file: Data cube of IRIS with its header (.fits file).

    Output:
    - Time vector as an array with elements in datetime.strptime format.

    This function allows us to create the time vector in datetime.strptime format for 
    IRIS data cubes. It takes the cadence from the header and calculates the time in 
    seconds for each frame in the data cube.The time vector is then created by adding 
    the calculated time to the start time of observation obtained from the header.
    The resulting time vector is returned as an array with elements in datetime.strptime format.
    """
    
    
    date_format   = "%Y-%m-%dT%H:%M:%S.%f"
    cadence       = iris_file[0].header['CADPL_AV']
    vec_sec_iris  = np.linspace(0,len(iris_file[0].data)-1,len(iris_file[0].data))*cadence 
    time_vec_iris = []
    for time_sec in vec_sec_iris:
        delta = timedelta(
        days=0,
        seconds=time_sec,
        microseconds=0,
        milliseconds=0,
        minutes=0,
        hours=0,
        weeks=0)
        time_vec_iris.append(delta+datetime.strptime(iris_file[0].header['STARTOBS'], date_format))
    return(np.array(time_vec_iris))

def fun_alma_rotate(cubo_alma, tilt_angle):
    """
    Remove the rotation from ALMA images and apply the aperture_circ_ALMA 
    function to match ALMA's FoV.

    Input:
    - cubo_alma: ALMA cube taken from SALSA.
    - tilt_angle: Tilted ALMA images.

    Output:
    - ALMA_CUBE without rotation.

    This function removes the rotation from ALMA images by rotating each 
    image in the cube using the 'SOLAR_P' angle from the header.
    The function then applies the aperture_circ_ALMA function to match 
    the Field of View (FoV) of ALMA.
    The resulting cube without rotation is returned.
    """
    #rotating ALMDA DATA
    imag_alma_rot=[]
    print('Rotating ALMA cube')
    for i in tqdm(range(len(cubo_alma)),colour='green'):
        imag_alma_rot.append(scnd.interpolation.rotate(cubo_alma[i],
                                                       angle=tilt_angle, reshape=False))
    print('rotate angle', tilt_angle)
    aperture=aperture_circ_ALMA(cubo_alma)
    imag_alma_rot                     = imag_alma_rot*aperture
    imag_alma_rot[imag_alma_rot<0.05] = cubo_alma[0][0][0] # cubo_alma[0][0][0] has mean value from ALMA cube
    return imag_alma_rot


def aperture_circ_ALMA(alma_cube):
    """
    Create a circular mask for a matrix with the Field of View (FoV) of ALMA.

    Input:
    - alma_cube: ALMA cube obtained using SALSA.

    Output:
    - A matrix mask with the same shape as ALMA's images.

    This function creates a circular mask for a matrix with the Field of View (FoV) of ALMA.
    The ALMA cube is used to determine the size of the mask.
    The mask is created by generating a grid of x and y values, calculating the radius of each point,
    and setting the values inside the circle to 1.0 in the aperture matrix.
    The resulting aperture matrix is returned.
    """
    x_size=alma_cube.shape[1]
    n = x_size
    #total size of the grid  in meters
    a = x_size
    aperture= np.zeros([n,n]) 
    #x grid 1D
    x = np.linspace(-a/2.0,a/2.0,n) 
    #y grid 1D
    y = np.linspace(-a/2.0,a/2.0,n)
    #radius of the telescope
    R = x_size/2 
    #making  circular aperture
    for j in range(0,n,1):
        for i in range(0,n,1):
            if ((x[i]**2.0 + y[j]**2.0)**(1/2) <= R):
                aperture[i,j] = 1.0
    return aperture

def find_nearest(array, alma_t0, alma_tf, name_telesc_of_array, name_telesc_t0_t_f ):
    '''
    Find the initial and final times in a time vector closest to the given ALMA observation times.

    Input:
    - array: Array with elements in datetime.datetime format (time vector obtained using the 
      DATA-OBS of each previously downloaded SDO frame).
    - alma_t0: Start time of the ALMA observation in datetime.datetime format.
    - alma_tf: End time of the ALMA observation in datetime.datetime format.
    - name_telesc_of_array: Name of the telescope corresponding to the array.
    - name_telesc_t0_t_f: Name of the telescope corresponding to alma_t0 and alma_tf.

    Output:
    - Positions in the vector array (array) closest to alma_t0 and alma_tf.

    This function allows us to find, in a time vector, the initial and final times at times alma_t0 
    and alma_tf.
    It calculates the indices of the elements in the array that are closest to alma_t0 and alma_tf.
    The function also prints the found times and their corresponding positions in the array.
    The indices idx_0 and idx_f are returned.
    '''
    array = np.asarray(array)
    idx_0 = (np.abs(array - alma_t0)).argmin()
    idx_f = (np.abs(array - alma_tf)).argmin()
    print('______________________________________________________')
    print('t_0_{} ='.format(name_telesc_of_array), array[idx_0])
    print('t_0_{}='.format(name_telesc_t0_t_f)   , alma_t0)
    print('______________________________________________________')
    print('t_f_{} ='.format(name_telesc_of_array), array[idx_f])
    print('t_0_{}='.format(name_telesc_t0_t_f)   , alma_tf)
    print('______________________________________________________')
    print('posicion_n_0:', idx_0, '| posicion_n_f:', idx_f)
    print('______________________________________________________')
    return idx_0, idx_f

def create_data_cube( path_img_crop, intrument):
    """
    This function create a date cube SDO AIA

    Last modified: Ordoñez Araujo, FJ  19.08.2022
    ------inputs----
    path_img     : it's the path where is the SDO image
    path_img_crop: it's the path where is the cropped SDO image
    
    intrument: this variable could be 'AIA' or 'HMI'
    and is necesary to differences on headers 
    -----outputs----
    data cube where its components is:
    
    data_cubo    :  data cube builed with SDO image  cropped 
    data_time_utc:  data_OBS of whole image  in forma %Y-%m-%dT%H:%M:%S.%f
    
    
    -----outputs----
    data_cube: Data cube built with the cropped SDO images.
    data_time_utc: OBS data of the entire image in the format %Y-%m-%dT%H:%M:%S.%f

    _______________________________________________________________________________________________________________
    WARNING: All matrices in the vector (data cube) should have the same shape.
    If the arrays in the data cube do not have the same shape, you will encounter problems when converting the list to an array.
    The error message "ValueError: could not broadcast input array from shape (453, 452) into shape (453,)" can occur.
    Please refer to the following URL for more information: 
    https://stackoverflow.com/questions/43977463/valueerror-could-not-broadcast-input-array-from-shape-224-224-3-into-shape-2
    ______________________________________________________________________________________________________________
    """
    intrument=intrument.lower()

    cubo_image_crop  = [] #cube of image crop
    date_OBS         = [] # data obs
    
    x_image_crop     = sorted(glob.glob(path_img_crop)) # of SDO image cropped

    #-----------------Creating a date cube------------------------
    """
    It is necessary to consider two cases because the DATE-OBS in HMI has 
    a different format ('2018-04-12T15:45:31.700001Z') than AIA 
    DATE-OBS ('2018-04-12T15:45:28.84') 
    """
    date_format = "%Y-%m-%dT%H:%M:%S.%f"
    
    size_x=[]
    size_y=[]
    if intrument=='hmi':
        print('|------Doing HMI cubic------|')
        print(len(x_image_crop))
        for i in tqdm(range(0,len(x_image_crop),1)):
            x=fits.open(x_image_crop[i])[0].data
            cubo_image_crop.append(x)
            
            time_1=fits.open(x_image_crop[i])[0].header['DATE-OBS'][0:-1]
            date_OBS.append(datetime.strptime(time_1, date_format))
            size_x.append(x.shape[1])
            size_y.append(x.shape[0])

    if intrument=='aia':
        print('|------Doing AIA cubic------|')
        print(len(x_image_crop))
        for i in tqdm(range(0,len(x_image_crop),1),desc='dirs'):
            x=fits.open(x_image_crop[i])[0].data
            cubo_image_crop.append(x) 
            
            time_1 =fits.open(x_image_crop[i])[0].header['DATE-OBS']
            date_OBS.append(datetime.strptime(time_1, date_format))
            size_x.append(x.shape[1])
            size_y.append(x.shape[0])
    
    #return np.array(cubo_image_crop), np.array(date_OBS) 
    if np.max(np.gradient(size_y))!= np.min(np.gradient(size_y)) or np.max(np.gradient(size_x))!= np.min(np.gradient(size_x)):
        print('----All snapshot does not has the same shape, we will make a small cut---')
        new_crop=[]
        for i in tqdm(range(0,len(cubo_image_crop),1)):
            new_crop.append(cubo_image_crop[i][0:int(np.min(size_y)),0:int(np.min(size_x))])
        return np.array(new_crop), np.array(date_OBS)  

    else:
        print('All frame has the same shape')
        return np.array(cubo_image_crop), np.array(date_OBS)  






def get_query_sdo(item, bottom_left, top_right, start_time, ends_time, email):
    """
    Generates and returns a search query object for solar data based on input parameters.
    The function evaluates whether the specified 'item' belongs to specific solar observation cadences
    or matches solar imaging instruments like 'HMI' or 'AIA'.

    The search criteria are determined by checking against pre-defined lists of cadences for the AIA instrument
    and identifying if the 'item' corresponds to HMI instruments.

    Parameters:
    - item (str): Identifier for the solar observation item. This can be an integer representing AIA wavelength 
                  (in angstroms) or a string ('hmi' or 'dopplergram') for HMI instruments.
    - bottom_left (tuple): Tuple of (x, y) coordinates representing the bottom-left corner of the desired 
                           helioprojective cutout region on the sun.
    - top_right (tuple): Tuple of (x, y) coordinates representing the top-right corner of the desired 
                         helioprojective cutout region.
    - start_time (str): The start time for the data query as a datetime object. 
    - ends_time (str): The end time for the data query as a datetime object.     
    # Time format should be ajusted such as, t_0 = str(2018-04-12 15:50:01), start_time = Time(t_0.isoformat(), scale='utc', format='isot')
    # Using from astropy.time import Time
    - email (str): Registered email address in the JSOC database which will be notified upon query completion.

    Returns:
    - query (QueryResponse object or None): A search query object configured for JSOC data retrieval via Fido. 
                                            Returns None if the 'item' does not match any pre-defined criteria.
	
    The query is built based on the specified 'item'. If the 'item' is found in either the aia_cad_12 or aia_cad_24 
    lists, the function configures a query with appropriate cadence, wavelength, and observational series. If 'item' 
    is 'hmi' or 'dopplergram', it adjusts the series and cadence accordingly. If no valid item is found, the function 
    returns None and prompts to check the provided information.
    """
    if item not in [94, 131, 171, 193, 211, 304, 335, 1600, 1700,  'hmi', 'aia', 'HMI', 'AIA', 'dopplergram', 'DOPPLERGRAM', 'Dopplergram']:
        raise ValueError("Supported methods are 'pearson' and 'ssim' only")
    # Define the cadence lists for AIA
    aia_cad_12 = [94, 131, 171, 193, 211, 304, 335]
    aia_cad_24 = [1600, 1700]

    # Check if the item is in the aia_cad_12 list
    if item in aia_cad_12:
        cadence = 12
        cutout = a.jsoc.Cutout(bottom_left, top_right=top_right, tracking=True)
        query = Fido.search(
            a.Time(start_time, ends_time),
            a.Wavelength(item*u.angstrom),
            a.Sample(cadence*u.second),
            a.jsoc.Series.aia_lev1_euv_12s,
            a.jsoc.Notify(email),
            a.jsoc.Segment.image,
            cutout)

    # Check if the item is in the aia_cad_24 list
    elif item in aia_cad_24:
        cadence = 24
        cutout = a.jsoc.Cutout(bottom_left, top_right=top_right, tracking=True)
        query = Fido.search(
            a.Time(start_time, ends_time),
            a.Wavelength(item*u.angstrom),
            a.Sample(cadence*u.second),
            a.jsoc.Series('aia.lev1_uv_24s'),
            a.jsoc.Notify(email),
            a.jsoc.Segment.image,
            cutout)

    # Check if the item matches 'hmi' (case insensitive)
    elif item.lower() == 'hmi':
        cutout = a.jsoc.Cutout(bottom_left, top_right=top_right, tracking=True)
        cadence = 45
        print('magnetogram, cadence:',cadence)
        query = Fido.search(
            a.Time(start_time, ends_time),
            a.Sample(cadence*u.second),
            a.jsoc.Series('hmi.M_45s'),
            a.jsoc.Notify(email),
            a.jsoc.Segment.magnetogram,
            cutout)
        
    elif item.lower() == 'dopplergram':
        cutout = a.jsoc.Cutout(bottom_left, top_right=top_right, tracking=True)
        cadence = 45
        print('dopplergram, cadence:',cadence)
        query = Fido.search(
            a.Time(start_time, ends_time),
            a.Sample(cadence*u.second),
            a.jsoc.Series('hmi.v_45s'),
            a.jsoc.Notify(email),
            a.jsoc.Segment.dopplergram,
            cutout)
    else:
        # Handle the case where the item is not in any list
        query = None
        print('Please check the provided information')
    return query


#---------------------------------------------------------------------------------------------------        
#---------------------------------------------------------------------------------------------------        
#---------------------------------------Alignment Functions ----------------------------------------
#---------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------
#---------------------------------- ALMA- IRIS -----------------------------------------------------
#---------------------------------------------------------------------------------------------------

def plot_alignment_results_alma_iris(mean_iris, mean_sdo, ccrr_pearson, x_max_pearson, y_max_pearson , radious,vmin_iris=None, vmax_iris=None):
    """
    Plots the alignment results including the average images from ALMA and IRIS and the correlation matrix.

    Parameters:
    - mean_iris (array): Average data from ALMA.
    - mean_sdo (array): Average data from IRIS.
    - ccrr_pearson (array): Pearson correlation matrix.
    - x_max_pearson (float): X-coordinate of the maximum Pearson correlation.
    - y_max_pearson (float): Y-coordinate of the maximum Pearson correlation.
    - radius (float): Radius for the circle to highlight max correlation.
    - vmin_iris (float, optional): Minimum data range for IRIS image plotting.
    - vmax_iris (float, optional): Maximum data range for IRIS image plotting.
    """
    fig, ax = plt.subplots(figsize=(13, 8), dpi=100, ncols=3, nrows=1, constrained_layout=True, gridspec_kw={'wspace': 0, 'hspace': 0})
    
    # Plotting the ALMA average data
    ax[0].imshow(mean_iris, cmap='hot', origin='lower')
    ax[0].set_title('ALMA Average')

    # Plotting the IRIS average data
    ax[1].imshow(mean_sdo, cmap='irissji2796', origin='lower', vmin=vmin_iris, vmax=vmax_iris)
    circle = matplotlib.patches.Circle((x_max_pearson, y_max_pearson), radious, ec='blue', fc='none', lw=1.5)
    ax[1].plot(x_max_pearson, y_max_pearson, 'b*', label='Maximum Correlation')
    ax[1].add_patch(circle)
    ax[1].set_title('IRIS Average')
    ax[1].legend()

    # Plotting the Pearson correlation matrix
    ax[2].imshow(ccrr_pearson, cmap='Greens', origin='lower')
    ax[2].plot(x_max_pearson, y_max_pearson, 'r*', label='Maximum Correlation')
    ax[2].set_title('Correlation Matrix')
    ax[2].legend()

    plt.show()

def calculate_correlation(iris_interp, alma_interp, alma_cube, method='pearson', scale_factor= 0.5545):
    """
    Calculates the correlation between interpolated IRIS and ALMA data using either Pearson or SSIM methods.

    Parameters:
    - iris_interp (array): Interpolated IRIS data.
    - alma_interp (array): Interpolated ALMA data.
    - alma_cube (array): ALMA data cube.
    - method (str): Method of correlation ('pearson' or 'ssim').
    - scale_factor (float): Scale factor for alignment. (IRIS PIXEL SIZE/ ALMA PIXEL SIZE)

    Returns:
    - tuple: Correlation matrix, max Y coordinate, max X coordinate, maximum correlation value.
    """
    method = method.lower()
    if method not in ['pearson', 'ssim']:
        raise ValueError("Supported methods are 'pearson' and 'ssim' only")

    # Initialize correlation matrices
    ccrr_pearson = np.zeros((iris_interp.shape[0], iris_interp.shape[1]))
    ccrr_ssim = np.zeros_like(ccrr_pearson)
    flatten_alma = alma_interp.flatten()

    if (iris_interp.shape[0] - alma_cube.shape[1])>160:
        counts_number = (iris_interp.shape[0] - alma_cube.shape[1])*(iris_interp.shape[1] - alma_cube.shape[2] - 2)
        print( f"Warning!  Perhaps SDO FoV is too large. You are calculating {counts_number} {method} coefficients.")

    print(f'|----Estimating {method.capitalize()} coefficients----|')
    # Calculate Pearson correlation or SSIM
    for i in tqdm(range(0, iris_interp.shape[0] - alma_cube.shape[1] - 2, 1), position=0, desc="Row", leave=False, colour='green', ncols=100):
        for j in range(0, iris_interp.shape[1] - alma_cube.shape[2] - 2, 1):
            matrix = iris_interp[i:i + alma_cube.shape[1], j:j + alma_cube.shape[2]].copy()
            x = i + alma_interp.shape[1] // 2
            y = j + alma_interp.shape[0] // 2

            if method == 'pearson':
                corr, _ = pearsonr(matrix.flatten(), flatten_alma)
                ccrr_pearson[y, x] = corr
            elif method == 'ssim':
                ssim = compare_ssim(matrix, alma_interp, win_size=11, data_range=matrix.max() - matrix.min())
                ccrr_ssim[y, x] = ssim
    
    # Find maximum correlation or SSIM location and scale them
    if method == 'pearson':
        result_matrix = ccrr_pearson
    else:
        result_matrix = ccrr_ssim

    y_max, x_max = np.unravel_index(np.argmax(result_matrix), result_matrix.shape)
    scaled_x_max = x_max / scale_factor
    scaled_y_max = y_max / scale_factor

    print('________________________________________________________________________')
    print(f'x_max_{method}={scaled_x_max} | y_max_{method}={scaled_y_max} | r={result_matrix.max()}')
    

    return result_matrix, scaled_y_max, scaled_x_max, result_matrix.max()


def align_iris_with_alma(iris_fits_file, timeutc_alma, alma_cube, alma_resolution, method='pearson', plot_results=True):
    """
    Aligns IRIS observations with ALMA data, scaling IRIS images to match ALMA resolution and computing correlations.

    Parameters:
    - iris_fits_file (fits object): FITS file containing IRIS data.
    - timeutc_alma (list): Time range for ALMA observations.
    - alma_cube (array): ALMA data cube.
    - alma_resolution (float): Resolution of ALMA data (arcsec).
    - method (str): Method to use for correlation analysis ('pearson' or 'ssim').
    - plot_results (bool): Whether to plot the results of the alignment.

    Returns:
    - tuple: Contains max X correlation, max Y correlation, initial position in IRIS array, final position in IRIS array, correlation matrix.
    """
    twave1 = iris_fits_file[0].header['TWAVE1']
    print(f'|-----------Aligning ALMA with IRIS {twave1} images----------|')
    iris_resolution = iris_fits_file[0].header['CDELT1']
    scale_factor = iris_resolution / alma_resolution  # to scale IRIS resolution
    print('scale_factor:', scale_factor)

    # Rescale IRIS image to match ALMA resolution
    print('Rescaling IRIS Image...')
    rescaling_iris = [rescale(iris_fits_file[0].data[i], scale_factor, anti_aliasing=False, order=0) for i in tqdm(range(len(iris_fits_file[0].data)), colour='green')]
    rescaling_iris = np.array(rescaling_iris)

    print('Old IRIS VECTOR', iris_fits_file[0].data.shape, 'New IRIS VECTOR', rescaling_iris.shape)

    # Creating time vector for IRIS and finding nearest matching time in ALMA data
    time_vec_iris = creating_time_vector(iris_fits_file)
    posicion_n_0, posicion_n_f = find_nearest(time_vec_iris, timeutc_alma[0], timeutc_alma[-1], 'IRIS', 'ALMA')

    print('________________________________________________________________________')
    print('t_0 position in the IRIS array:', int(posicion_n_0), '| t_f position in the IRIS array:', int(posicion_n_f))

    # Compute mean data for IRIS and ALMA within the selected time window
    mean_iris = np.mean(rescaling_iris[int(posicion_n_0):int(posicion_n_f)], axis=0)
    mean_alma = np.mean(alma_cube, axis=0)
    iris_interp = np.interp(mean_iris, (np.percentile(mean_iris, 12), np.percentile(mean_iris, 95)), (0, 1))
    alma_interp = np.interp(mean_alma, (np.percentile(mean_alma, 11), np.percentile(mean_alma, 95)), (0, 1))

    vmin_iris = np.percentile(mean_iris, 11)
    vmax_iris = np.percentile(mean_iris, 95)

    # Perform correlation analysis
    radius = mean_alma.shape[1] / 2
    correlation_matrix, y_max_correlation, x_max_correlation, _ = calculate_correlation(iris_interp, alma_interp, alma_cube, method=method, scale_factor=scale_factor)

    # Optionally plot the results
    if plot_results==True:
        try:
            plot_alignment_results_alma_iris(mean_alma, mean_iris, correlation_matrix, x_max_correlation * scale_factor, y_max_correlation * scale_factor, radius, vmin_iris, vmax_iris)
        except :
            print('There was an error, we could not plot results:')

    return x_max_correlation, y_max_correlation, posicion_n_0, posicion_n_f, correlation_matrix




#---------------------------------------------------------------------------------------------------
#---------------------------------- IRIS- SDO ------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def plot_alignment_results_iris_sdo(mean_iris, mean_sdo, ccrr_pearson, x_max_pearson, y_max_pearson ,vmin_iris=None, vmax_iris=None):
    """
    Plots the alignment results including the average images from ALMA and IRIS and the correlation matrix.

    Parameters:
    - mean_iris (array): Average data from ALMA.
    - mean_sdo (array): Average data from IRIS.
    - ccrr_pearson (array): Pearson correlation matrix.
    - x_max_pearson (float): X-coordinate of the maximum Pearson correlation.
    - y_max_pearson (float): Y-coordinate of the maximum Pearson correlation.
    - radius (float): Radius for the circle to highlight max correlation.
    - vmin_iris (float, optional): Minimum data range for IRIS image plotting.
    - vmax_iris (float, optional): Maximum data range for IRIS image plotting.
    """
    vmin_iris = np.percentile(mean_iris, 11)
    vmax_iris = np.percentile(mean_iris, 95)


    fig, ax = plt.subplots(figsize=(13, 8), dpi=100, ncols=3, nrows=1, constrained_layout=True, gridspec_kw={'wspace': 0, 'hspace': 0})
    
    # Plotting the ALMA average data
    ax[0].imshow(mean_iris, cmap='irissji2796', origin='lower', vmin=vmin_iris, vmax=vmax_iris)
    ax[0].set_title('IRIS Average')

    # Plotting the IRIS average data
    ax[1].imshow(mean_sdo, cmap='irissji1330', origin='lower')
    ax[1].plot(x_max_pearson, y_max_pearson, 'b*', label='Maximum Correlation')
    ax[1].set_title('SDO Average')
    ax[1].legend()

    # Plotting the Pearson correlation matrix
    ax[2].imshow(ccrr_pearson, cmap='Greens', origin='lower')
    ax[2].plot(x_max_pearson, y_max_pearson, 'r*', label='Maximum Correlation')
    ax[2].set_title('Correlation Matrix')
    ax[2].legend()

    plt.show()

def plot_alignment_results_alma_sdo(mean_iris, mean_sdo, ccrr_pearson, x_max_pearson, y_max_pearson ,vmin_iris=None, vmax_iris=None):
    """
    Plots the alignment results including the average images from ALMA and IRIS and the correlation matrix.

    Parameters:
    - mean_iris (array): Average data from ALMA.
    - mean_sdo (array): Average data from IRIS.
    - ccrr_pearson (array): Pearson correlation matrix.
    - x_max_pearson (float): X-coordinate of the maximum Pearson correlation.
    - y_max_pearson (float): Y-coordinate of the maximum Pearson correlation.
    - radius (float): Radius for the circle to highlight max correlation.
    - vmin_iris (float, optional): Minimum data range for IRIS image plotting.
    - vmax_iris (float, optional): Maximum data range for IRIS image plotting.
    """
    vmin_iris = np.percentile(mean_iris, 11)
    vmax_iris = np.percentile(mean_iris, 95)


    fig, ax = plt.subplots(figsize=(13, 8), dpi=100, ncols=3, nrows=1, constrained_layout=True, gridspec_kw={'wspace': 0, 'hspace': 0})
    
    # Plotting the ALMA average data
    ax[0].imshow(mean_iris, cmap='hot', origin='lower')
    ax[0].set_title('ALMA Average')

    # Plotting the IRIS average data
    ax[1].imshow(mean_sdo, cmap='irissji1330', origin='lower')
    ax[1].plot(x_max_pearson, y_max_pearson, 'b*', label='Maximum Correlation')
    ax[1].set_title('SDO Average')
    ax[1].legend()

    # Plotting the Pearson correlation matrix
    ax[2].imshow(ccrr_pearson, cmap='Greens', origin='lower')
    ax[2].plot(x_max_pearson, y_max_pearson, 'r*', label='Maximum Correlation')
    ax[2].set_title('Correlation Matrix')
    ax[2].legend()

    plt.show()




def calculate_correlation_pearson(sdo_img, iris_img):
    """
    Calculate the Pearson correlation coefficient matrix for two images.
    
    Parameters:
        sdo_img (np.ndarray): The SDO image array.
        iris_img (np.ndarray): The IRIS image array resized to SDO scale.
    
    Returns:
        np.ndarray: Matrix of Pearson correlation coefficients.
    """
    if (sdo_img.shape[0] - iris_img.shape[0]+1)>150:
        counts_number = (sdo_img.shape[0] - iris_img.shape[0]+1)*(sdo_img.shape[1] - iris_img.shape[1]+1)
        print( f"Warning! Perhaps SDO FoV is too large. You are calculating {counts_number} Pearson coefficients.")
    correlation_matrix = np.zeros((sdo_img.shape[0], sdo_img.shape[1]))
    for i in tqdm(range(0, sdo_img.shape[0] - iris_img.shape[0]+1, 1), position=0, desc="i", leave=False, colour='green', ncols=100):
        for j in  range(0, sdo_img.shape[1] - iris_img.shape[1]+1, 1):
            sub_sdo_img = sdo_img[0+i:iris_img.shape[0]+i, 0+j:iris_img.shape[1]+j].copy()
            corr, _ = pearsonr(sub_sdo_img.flatten(), iris_img.flatten())

            x = i+iris_img.shape[1]//2
            y = j+iris_img.shape[0]//2
            correlation_matrix[x, y] = corr
    return correlation_matrix

def calculate_ssim(sdo_img, iris_img):
    """
    Calculate the Structural Similarity Index (SSIM) for each position where the IRIS image
    can fully overlay the SDO image.
    
    Parameters:
        sdo_img (np.ndarray): The SDO image array.
        iris_img (np.ndarray): The IRIS image array resized to SDO scale.
    
    Returns:
        np.ndarray: Matrix of SSIM values.
    """
    if (sdo_img.shape[0] - iris_img.shape[0]+1)>150:
        counts_number = (sdo_img.shape[0] - iris_img.shape[0]+1)*(sdo_img.shape[1] - iris_img.shape[1]+1)
        print( f"Warning! Perhaps SDO FoV is too large. You are calculating {counts_number} SSIM coefficients.")

    ssim_matrix = np.zeros((sdo_img.shape[0], sdo_img.shape[1]))
    for i in tqdm(range(0, sdo_img.shape[0] - iris_img.shape[0]+1, 1), position=0, desc="i", leave=False, colour='green', ncols=100):
        for j in  range(0, sdo_img.shape[1] - iris_img.shape[1]+1, 1):
            sub_sdo_img = sdo_img[0+i:iris_img.shape[0]+i, 0+j:iris_img.shape[1]+j].copy()
            corr, _ = pearsonr(sub_sdo_img.flatten(), iris_img.flatten())

            x = i+iris_img.shape[1]//2
            y = j+iris_img.shape[0]//2
            ssim_matrix[x, y] = corr
    return ssim_matrix

def align_iris_with_sdo(sdo_file, iris_data, iris_resolution, time_vec_iris, time_vec_sdo, time_vec_alma, path_img_crop, method='pearson', plot_results=True):
    """
    Aligns IRIS image data with SDO image data by rescaling the IRIS images and calculating
    the correlation or similarity with corresponding SDO images to find the best alignment.

    Parameters:
        sdo_file (np.ndarray): Array of frames from the SDO/AIA, loaded from a .npy file.
        iris_data (np.ndarray): Data cube from the IRIS database.
        iris_resolution (float): Spatial resolution of the IRIS images.
        time_vec_iris (np.ndarray): Timestamps for each frame in the IRIS data cube.
        time_vec_sdo (np.ndarray): Timestamps for each frame in the SDO data cube.
        time_vec_alma (np.ndarray): Timestamps for each frame in the ALMA data cube.
        path_img_crop (str): File path to the directory containing cropped images.
        method (str): Method used for alignment ('pearson', 'ssim'). Default is 'pearson'.
        plot_results (bool): If True, plot the alignment results. Default is True.

    Returns:
        tuple: Estimated helioprojective coordinates (Tx, Ty) of the center of the IRIS matrix.
    """
    method = method.lower()
    if method not in ['pearson', 'ssim']:
        raise ValueError("Supported methods are 'pearson' and 'ssim' only")

    image_files = sorted(glob.glob(path_img_crop))
    sdo_resolution = fits.open(image_files[0])[0].header['CDELT2']
    scale_factor = iris_resolution / sdo_resolution
    print(f'Rescaling IRIS images with scale factor: {scale_factor}')

    # Rescale IRIS images
    rescaled_iris = [rescale(img, scale_factor, order=0, anti_aliasing=False, preserve_range=True) for img in iris_data]
    rescaled_iris = np.array(rescaled_iris)

    print('Initial shape:', sdo_file.shape, '| Rescaled shape:', rescaled_iris.shape)

    # Align the time vectors with ALMA
    start_idx_sdo, end_idx_sdo   = find_nearest(time_vec_sdo  , time_vec_alma[0], time_vec_alma[-1], 'SDO', 'IRIS')
    start_idx_iris, end_idx_iris = find_nearest(time_vec_iris, time_vec_alma[0], time_vec_alma[-1], 'SDO', 'IRIS')

    # Compute mean images in the aligned time window
    mean_sdo = np.mean(sdo_file[start_idx_sdo:end_idx_sdo], axis=0)
    mean_iris = np.mean(rescaled_iris[start_idx_iris:end_idx_iris], axis=0)

    # Normalize images
    sdo_normalized = np.interp(mean_sdo, (np.percentile(mean_sdo, 10), np.percentile(mean_sdo, 90)), (0, 1))
    iris_normalized = np.interp(mean_iris, (np.percentile(mean_iris, 10), np.percentile(mean_iris, 90)), (0, 1))

    # Calculate alignment using the specified method
    if method.lower() == 'pearson':
        correlation_matrix = calculate_correlation_pearson(sdo_normalized, iris_normalized)
        y_max, x_max = np.unravel_index(np.argmax(correlation_matrix), correlation_matrix.shape)
    elif method.lower() == 'ssim':
        similarity_matrix = calculate_ssim(sdo_normalized, iris_normalized)
        y_max, x_max = np.unravel_index(np.argmax(similarity_matrix), similarity_matrix.shape)

    # Convert pixel maxima to helioprojective coordinates
    map_sunpy = sunpy.map.Map(image_files[start_idx_sdo])
    helioprojective_coords = map_sunpy.pixel_to_world(x_max * u.pix, y_max * u.pix)
    print(f'Best alignment (Tx, Ty) in arcsec: {helioprojective_coords.Tx.value}, {helioprojective_coords.Ty.value}')

    # Optionally plot results
    if plot_results == True:
        try:
            plot_alignment_results_iris_sdo(mean_iris, mean_sdo, correlation_matrix if method == 'pearson' else similarity_matrix, x_max, y_max)
        except :
            print('There was an error, we could not plot results:')

    return helioprojective_coords.Tx.value, helioprojective_coords.Ty.value



#---------------------------------------------------------------------------------------------------
#---------------------------------- ALMA- SDO ------------------------------------------------------
#---------------------------------------------------------------------------------------------------

def align_alma_with_sdo(sdo_file, time_filter, time_alma, alma_cube, alma_resolution, path_mapMap, method='pearson', plot_results=True):
    """
    Aligns Solar Dynamics Observatory (SDO) data with Atacama Large Millimeter/submillimeter Array (ALMA) observations.

    This function aligns data from SDO and ALMA using cross-correlation techniques specified by the method argument. 
    It rescales the resolution of SDO data to match that of ALMA, selects the relevant observation window, 
    and computes the alignment using either Pearson correlation or Structural Similarity Index (SSIM).

    Parameters:
        sdo_file (list or array): Collection of SDO images to be processed.
        time_filter (list)      : Two-element list defining the start and end times for filtering the observations.
        time_alma (list)        : Time stamps corresponding to ALMA observations to align with the SDO data.
        alma_cube (np.array)    : Data cube containing ALMA observations.
        alma_resolution (float) : Spatial resolution of the ALMA observations.
        path_mapMap (str)       : Path pattern to locate FITS files containing SDO data.
        method (str, optional)  : Method to use for calculating alignment; supported values are 'pearson' for Pearson 
        correlation, or 'ssim' for Structural Similarity Index. Default is 'pearson'.

    Returns:
        tuple: Contains helioprojective coordinates (Tx, Ty in arcseconds) of the alignment 
        point, and the starting and ending indexes of the selected observation window in SDO data.
    """

    method = method.lower()
    if method not in ['pearson', 'ssim']:
        raise ValueError("Supported methods are 'pearson' and 'ssim' only")

    url = sorted(glob.glob(path_mapMap))

  
    
    try:
        name_filter = fits.open(url[0])[0].header['TELESCOP']+ \
                  str('-')+fits.open(url[0])[0].header['WAVE_STR']
        print('|-----------------{}-------------------|'.format(name_filter))  
      
    except:
        name_filter = fits.open(url[0])[0].header['TELESCOP']+ \
              str('-')+fits.open(url[0])[0].header['INSTRUME']
        print('|-----------------{}-------------------|'.format(name_filter)) 
 
    # To scale IRIS resolution
    sdo_resolution = fits.open(url[0])[0].header['CDELT1']
    cadence_filter = time_filter[1]-time_filter[0]
    scale_factor = sdo_resolution/alma_resolution
    print('scale_factor:', scale_factor)
    print('Rescaling AIA  images')
    # For example to AIA 94A filter
    # from (226, 152, 151) ----> (950, 429, 408)
    # same resolution to SDO/AIA AND HMI-ALMA 0.3
    rescaled_sdo = []
    for i in tqdm(range(len(sdo_file)), colour='green'):
        rescaled_sdo.append(rescale(sdo_file[i], scale_factor, order=0,
                                    anti_aliasing=False, preserve_range=True))
    rescaled_sdo = np.array(rescaled_sdo)

    print('________________________________________________________________________')
    print('initial shape:', sdo_file.shape,
          '| rescaled shape:', rescaled_sdo.shape)
    # --------------Selecting the same windown of ALMA----------------
    posicion_n_0, posicion_n_f = find_nearest( time_filter, time_alma[0], time_alma[-1], 'ALMA', 'SDO')

    # -------------- Mean IRIS and ALMA datacube-----------------------
    mean_sdo    = np.sum([i for i in rescaled_sdo[int(posicion_n_0):int(
                          posicion_n_f)]], 0)/len(rescaled_sdo[int(posicion_n_0):int(posicion_n_f)])

    mean_alma   = np.sum   ([i for i in alma_cube], 0)/len(alma_cube)
    sdo_interp  = np.interp(mean_sdo , (np.percentile(mean_sdo, 11) , np.percentile(mean_sdo , 95)), (0, +1))
    alma_interp = np.interp(mean_alma, (np.percentile(mean_alma, 11), np.percentile(mean_alma, 95)), (0, +1))

    # Calculate alignment using the specified method
    if method.lower() == 'pearson':
        correlation_matrix = calculate_correlation_pearson(sdo_interp, alma_interp)
        y_max, x_max = np.unravel_index(np.argmax(correlation_matrix), correlation_matrix.shape)
    elif method.lower() == 'ssim':
        similarity_matrix = calculate_ssim(sdo_interp, alma_interp)
        y_max, x_max = np.unravel_index(np.argmax(similarity_matrix), similarity_matrix.shape)
    

    if plot_results==True:
        try:
            plot_alignment_results_alma_sdo(mean_alma, mean_sdo, correlation_matrix, x_max, y_max, vmin_iris=None, vmax_iris=None)
        except :
            print('There was an error, we could not plot results:')


    print('________________________________________________________________________')
    print(f'x_max_{method}=', x_max, f'| y_max_{method}=',
          y_max, "| r=", np.max(correlation_matrix))

    print('________________________________________________________________________')
    print('|-----Helioprojective Coordinates-------|')
    map_sunpy                            = sunpy.map.Map(url[posicion_n_0])
    x_scaled , y_scaled  = x_max/scale_factor   ,  y_max/scale_factor
    coordinates                         = map_sunpy.pixel_to_world(x_scaled *u.pix,  y_scaled *u.pix )

    print('Pearson: (Tx, Ty) in arcsec',coordinates.Tx.value,coordinates.Ty.value)
    return coordinates.Tx.value, coordinates.Ty.value, posicion_n_0, posicion_n_f     



