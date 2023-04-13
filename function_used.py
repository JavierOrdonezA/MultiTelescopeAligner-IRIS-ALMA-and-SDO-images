from astropy.coordinates import SkyCoord
from skimage.transform import rescale
from sunpy.net import Fido, attrs as SAF
from  scipy.stats.stats import pearsonr
from skimage.metrics import structural_similarity as ssim
from datetime import  date, datetime,timedelta
import scipy.ndimage as scnd
from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np 
import sunpy.map
import numpy as np 
import sunpy
import glob 
import os
import warnings
warnings.filterwarnings("ignore")

def download_sdo(tsdo_init_dwl,tsdo_fina_dwl,Tx, Ty, aiasq):
    """
    This function download and crop SDO cubes

    Last modified: F.J Ordoñez A,  7.01.2023 

    Parameters:
    -----------
    tsdo_init_dwl: datetime.datetime
        Starting time 
    tsdo_fina_dwl: datetime.datetime
        Ennding time
    Tx: float 
        Centre position in X   helio projective coordinates 
    Ty: float 
        Centre position in Y   helio projective coordinates 
    aiasq: float
        widht of desired region
    Returns:
    ----------

    It stores the downloaded data in the corresponding folders 

    Comments:
    ---------
    It can be extended ot oher wavelenghts
    """
    
    ## AIA wavelength
    wl=[94,131,171,193,211,304,335,1600,1700]
   

    # Make Floders to store downloads.
    os.system('mkdir SDO')
    for i in range(len(wl)):
        os.system('mkdir SDO/{}'.format(wl[i]))
        os.system('mkdir SDO/{}/crop'.format(wl[i]))
        os.system('mkdir SDO/{}/crop_and_rotate'.format(wl[i]))      
    
    os.system('mkdir SDO/magt')
    os.system('mkdir SDO/magt/crop')
    os.system('mkdir SDO/magt/crop_and_rotate')
    
    ##Downloading AIA
    
    for item in tqdm(wl):
        print('Downloading %i A'%item)
        result = Fido.search(SAF.Time(tsdo_init_dwl, tsdo_fina_dwl), 
                                 SAF.Instrument("aia"), SAF.Wavelength(item*u.angstrom))
        print(result)
        file_download = Fido.fetch(result, path='SDO/{}'.format(item))
        
        if len(result)==0:
            continue
        
        if len(file_download)!=len(result[0]):
            
            print('There was a error in server downloading the data')
            
        while len(file_download)!=len(result[0]):
            
            print('Retrying the download, this is normal')
            
            file_download = Fido.fetch(result, path='SDO/{}'.format(item))
    
        #Saving AIA time sequence
        file_download =sorted(file_download)
        # Load all files into a Map sequence
        tmp = sunpy.map.Map(file_download)
        aia_seq = []
        # cropping 
        print('doing cropping')
        for img in tmp:
            top_right = SkyCoord((Tx+aiasq)*u.arcsec, (Ty+aiasq)*u.arcsec, frame=img.coordinate_frame)
            bottom_left = SkyCoord((Tx-aiasq)*u.arcsec, (Ty-aiasq)*u.arcsec, frame=img.coordinate_frame)
            aia_seq.append(img.submap(bottom_left, top_right=top_right)) 


        aia_seq_totated=[]
        for i in range(0,len(aia_seq),1):    
            img= aia_seq[i].rotate(use_scipy=True, missing=0)
            aia_seq_totated.append(img)


        new_files_totated = [str('SDO/{}/crop_and_rotate/').format(item)+file_download[a][7:-6]+str('_crop_and_rotate.fits') for a in range(len(file_download))]
    
        for img, file_name in zip(aia_seq_totated, new_files_totated):
            img.save(file_name)


        new_files = [str('SDO/{}/crop/').format(item)+file_download[a][7:-6]+str('_crop.fits') for a in range(len(file_download))]
        for img, file_name in zip(aia_seq, new_files):
            img.save(file_name)
            
    del wl,result, file_download, tmp, aia_seq, new_files, img, aia_seq_totated, new_files_totated
   #--------------------------------------------------------------------------------------------------
    
    ###Downloading Magnetograms
    wl = ['LOS_magnetic_field']
    print('Donwloading SDO/HMI data')
    for item in tqdm(wl):
        print('Downloading %s'%item)
        result = Fido.search(SAF.Time(tsdo_init_dwl, tsdo_fina_dwl), 
                             SAF.Instrument("hmi"), SAF.Physobs(item), 
                             SAF.Sample(45*u.second))    
        print(result)
        file_download = Fido.fetch(result, path='SDO/magt',site='ROB')
        
        if len(file_download)!=len(result[0]):    
            print('There was a error in server')
            
        while len(file_download)!=len(result[0]):
            print('Continue the download to alternative form')
            file_download = Fido.fetch(result, path='SDO/magt',site='ROB')
        
        
        #Saving HMI time sequence
        file_download = sorted(file_download)
        ###Load all files in map sequence
        tmp = sunpy.map.Map(file_download)
        #Cropping defined area
        
        # In this case, we do not run aiaprep because we're only using one channel
        aia_seq = []  #It could be called hmi_seq but it is just a variable
        print('cropping img')
        for img in tmp:
            top_right = SkyCoord((Tx+aiasq)*u.arcsec, (Ty+aiasq)*u.arcsec, frame=img.coordinate_frame)
            bottom_left = SkyCoord((Tx-aiasq)*u.arcsec, (Ty-aiasq)*u.arcsec, frame=img.coordinate_frame)
            aia_seq.append(img.submap(bottom_left, top_right=top_right))


        aia_seq_totated=[]
        for i in range(0,len(aia_seq),1):    
            img= aia_seq[i].rotate(use_scipy=True, missing=0)
            aia_seq_totated.append(img)

        new_files_totated = [str('SDO/magt/crop_and_rotate/')+file_download[a][8:-5]+str('_crop_and_rotate.fits') for a in range(len(file_download))]
        for img, file_name in zip(aia_seq_totated, new_files_totated):
            img.save(file_name)
        
        new_files =  [str('SDO/magt/crop/')+file_download[a][8:-5]+str('_crop.fits') for a in range(len(file_download))]
        for img, file_name in zip(aia_seq, new_files):
            img.save(file_name)
    return print('The download process for SDO and HMIhasve finished.')





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
    _______________________________________________________________________________________________________________
    WARNING:, ALL matrixs in the vector (data cube) should has the same shape, in general 
    when we do  the crop using sunpy this this will not be.
    when all the arrays in the data cube do not have the same shape you will have 
    problems trying to convert the list to an array . The error (likes, ValueError: could not 
    broadcast input array from shape (453,452) into shape (453,)) you will get was 
    reported in 
    url : https://stackoverflow.com/questions/43977463/valueerror-could-not-broadcast-input-array-from-shape-224-224-3-into-shape-2
    ______________________________________________________________________________________________________________
    """
    intrument=intrument.lower()

    cubo_image_crop  = [] #cube of image crop
    date_OBS         = [] # data obs
    
    x_image_crop     = sorted(glob.glob(path_img_crop)) # of SDO image cropped

    #-----------------Creating a date cube------------------------
    """
     It necesary considere two case becasuse the DATE-OBS in HMI has 
     another format ('2018-04-12T15:45:31.700001Z') that  AIA 
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


def aperture_circ_ALMA(alma_cube):
    """
    This function allow us to do a mask of a matrix with the FoV of ALMA

    Input  : almacube obtained using SALSA 
    Output : A matrix mask with the same shape to ALMA'S IMAGES 
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

def fun_alma_rotate(cubo_alma, header_alma):
    """
    Regard to IRIS/SDO images the ALMA images are rotated. 
    This function delete this rotation and applied  the 
    function aperture_circ_ALMA to has the FoV of ALMA
    
    Input : ALMA CUBE taken from  SALSA
    Output: ALMA_CUBE without rotation  
    """
    #rotating ALMDA DATA
    imag_alma_rot=[]
    print('Rotating ALMA cube')
    for i in tqdm(range(len(cubo_alma)),colour='green'):
        imag_alma_rot.append(scnd.interpolation.rotate(cubo_alma[i],
                                                       angle=header_alma['SOLAR_P'],reshape=False))
    print('rotate angle',header_alma['SOLAR_P'])
    aperture=aperture_circ_ALMA(cubo_alma)
    imag_alma_rot                     = imag_alma_rot*aperture
    imag_alma_rot[imag_alma_rot<0.05] = cubo_alma[0][0][0] # cubo_alma[0][0][0] has mean value from ALMA cube
    return imag_alma_rot






def iris_rescale(iris_fits_file, alma_resolution):
    """
    This function rescale IRISeS images. For that,  
    we  taken  the IRIS'S images in its resolution 
    and we carry  to the ALMA resolution

    Input  :  IRIS FILE.FIT 
    Output :  Vector with IRIS IMAGES  rescaled 
    """
    # wearing  IRIS' resolution to ALMA resolution
    #|---- Band 3 ALMA resolution is 0.3 arcsec----|
    #|------- header['CDELT1A']=0.3 arcsec---------|
    print('______________________________________________________')
    print('initial shape:',iris_fits_file[0].data.shape)

    #alma_resolution=0.3 #|---- Band 3 ALMA resolution is 0.3 arcsec----|
    iris_resolution=iris_fits_file[0].header['CDELT1']
    scale_factor=iris_resolution/alma_resolution # to scale IRIS resolution
    print('scale_factor:', scale_factor)
    rescaling_iris=[]

    # from (950, 774, 735) ----> (950, 429, 408)
    #same resolution to IRIS-ALMA 0.3
    for i  in  tqdm(range(len(iris_fits_file[0].data)),colour='green' ):
        rescaling_iris.append(rescale(iris_fits_file[0].data[i], scale_factor, anti_aliasing=False,order=0,preserve_range=True))
    rescaling_iris=np.array(rescaling_iris) 
    print('______________________________________________________')
    print('final shape:',rescaling_iris.shape)
    
    return rescaling_iris

def creating_time_vector(iris_file):
    """
    This function allows us to create the time vector 
    in datetime.strptime format for IRIS data cubes.
    ----------------------------Input----------------------------------
    Parameters: 
                iris_file: type .fits  is the datacube of iris with its header 

    ----------------------------Output----------------------------------
    time vector type array and its elements  is  datetime.strptime format 
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

def find_nearest(array, alma_t0, alma_tf, name_telesc_of_array, name_telesc_t0_t_f ):
    '''
    This function allows us to find (in a time vector) 
    the initial and final times at times alma_t0, alma_tf.
    
    --------------------Input----------------------------------
    Parameters: 
                array : type array and its elements is in the datetime.datetime 
                        format (in this case we  obtained this vector using the )
                        in our case this time vector was obtained with the DATA-OBS 
                        of each previously downloaded SDO frame  
                alma_t0: type T0 is  where the ALMA observation has started
                         datetime.datetime  format
                alma_tf: type Tf is  where the ALMA observation has finalized 
                         datetime.datetime  format  
    --------------------Output----------------------------------
    Positions in the vector array (array) closest to the   alma_t0, alma_tf        
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



def similarity(iris_fits_file, timeutc_alma, alma_cube, alma_resolution, method='pearson'):
    """
    This function aligns the IRIS and ALMA image. For that, it estimates the best Pearson 
    coefficient between the average matrix of the ALMA and IRIS data in the same observation 
    window and on the same scale (imsolating the images between a minimum (12th Percentile) 
    and a maximum (95th Percentile)). and for a given center, the ALMA image is transposed 
    onto the IRIS image in search of a maximum Pearson correlation.

    inpunt: -IRIS FILE
            -TIME UTC and ALMACUBE taken from  SALSA   
            -alma_resolution : pixel size of alma (in band 3 is 0.3 arcsec)
    Output: - x_max_pearson    : Position in the IRIS matrix that has the best x correlation 
            - y_max_pearson    : Position in the IRIS matrix that has the best y correlation 
            - int(posicion_n_0): Position in the IRIS datacube of image that was taken in the 
                                 t_0 ALMAUTC
            -int(posicion_n_f)   Position in the IRIS datacube of image that was taken in the 
                                 t_f ALMAUTC
    
    """
    print('|-----------Aligning ALMA  with IRIS {}  images----------| '.format(iris_fits_file[0].header['TWAVE1']))
    #| wearing  IRIS' resolution to ALMA resolution
    #|---- Band 3 ALMA resolution is 0.3 arcsec----|
    #|------- header['CDELT1A']=0.3 arcsec---------|

    #|---- Band 3 ALMA resolution is 0.3 arcsec----|
    iris_resolution=iris_fits_file[0].header['CDELT1']
    scale_factor=iris_resolution/alma_resolution # to scale IRIS resolution
    print('scale_factor:', scale_factor)
    rescaling_iris=[]

    # from (950, 774, 735) ----> (950, 429, 408)
    #same resolution to IRIS-ALMA 0.3
    print('ReScaling IRIS Image')
    for i  in tqdm(range(len(iris_fits_file[0].data)), colour='green'):
        rescaling_iris.append(rescale(iris_fits_file[0].data[i], scale_factor, anti_aliasing=False, order=0))
    rescaling_iris=np.array(rescaling_iris) 
    
    print(print('Old IRIS VECTOR', iris_fits_file[0].data.shape),'New IRIS VECTOR', rescaling_iris.shape)
    #print(np.max(rescaling_iris),np.min(rescaling_iris), np.std(rescaling_iris))
    #--------------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------------
    
    #--------------Selecting the same windown of ALMA----------------

    time_vec_iris = creating_time_vector(iris_fits_file)

    posicion_n_0 , posicion_n_f  = find_nearest(time_vec_iris, timeutc_alma[0], timeutc_alma[-1], 'IRIS', 'ALMA')


    print('________________________________________________________________________')
    print('t_0 position in the IRIS array:',int(posicion_n_0),'| t_f position in the IRIS array:',int(posicion_n_f))
    #-----------------------------------------------------------------
    
    #-------------- Mean IRIS and ALMA datacube-----------------------
    mean_iris_iris = np.sum([i for i in rescaling_iris[int(posicion_n_0):int(posicion_n_f)]],0)/len(rescaling_iris[int(posicion_n_0):int(posicion_n_f)])
    mean_alma      = np.sum([i for i in alma_cube] , 0)/len(alma_cube)
    iris_interp    = np.interp(mean_iris_iris , (np.percentile(mean_iris_iris,12), np.percentile(mean_iris_iris,95)), (0,+1))
    alma_interp    = np.interp(mean_alma      , (np.percentile(mean_alma,11)     , np.percentile(mean_alma,95))     , (0,+1))


    #-------------- PEARSON ANALYSIS ---------------------------------    
    ccrr_pearson = np.zeros((iris_interp.shape[0],iris_interp.shape[1]))
    ccrr_ssim = np.zeros((iris_interp.shape[0], iris_interp.shape[1]))
    
    # ALMA_VECTOR
    flatten_alma=alma_interp.flatten()
    method = method.lower()
    
        
    if method== 'pearson':
        print('|----Estimating Pearson coeffients---- | ')
        for i in tqdm(range(0,iris_interp.shape[0]-255,1), position=0, desc="i", leave=False, colour='green', ncols=100):
            for j in range(0,iris_interp.shape[1]-255,1):
                matrix=iris_interp[0+i:253+i,0+j:253+j].copy()
                #matrix=matrix.flatten()
                corr, valor_p= pearsonr(matrix.flatten(), flatten_alma)
          #      ssim_none = ssim(matrix, alma_interp,
          #                       data_range=iris_interp.max() - iris_interp.min())

                x=i+alma_interp.shape[1]//2
                y=j+alma_interp.shape[0]//2

                ccrr_pearson[x][y]=corr
         #       ccrr_ssim[x][y] = ssim_none
        
        x_max_pearson = np.where(ccrr_pearson==np.max(ccrr_pearson))[1][0]
        y_max_pearson = np.where(ccrr_pearson==np.max(ccrr_pearson))[0][0]
        
        print('________________________________________________________________________')
        print('x_max_pearson=',x_max_pearson  , '| y_max_pearson=', y_max_pearson,"| r=",np.max(ccrr_pearson))
        return  x_max_pearson/scale_factor, y_max_pearson/scale_factor, posicion_n_0 , posicion_n_f 
        #-----------------------------------------------------------------------------------------------------------
        
        if method== 'ssim':
            print('|----Estimating Structural Similarity Index  (SSIM)---- | ')
    
            print('WARNING: Estimating the similarity between the images using this method is computationally expensive.')
            for i in tqdm(range(0,iris_interp.shape[0]-255,1), position=0, desc="i", leave=False, colour='green', ncols=100):
                for j in range(0,iris_interp.shape[1]-255,1):
                    matrix=iris_interp[0+i:253+i,0+j:253+j].copy()
                    
                    
                    ssim_none = ssim(matrix, alma_interp,
                                     data_range=iris_interp.max() - iris_interp.min())

                    x=i+alma_interp.shape[1]//2
                    y=j+alma_interp.shape[0]//2

                    ccrr_ssim[x][y] = ssim_none
            


            x_max_ssim = np.where(ccrr_ssim == np.max(ccrr_ssim))[1][0]
            y_max_ssim = np.where(ccrr_ssim == np.max(ccrr_ssim))[0][0]
            print('________________________________________________________________________')
            print('x_max_ssim   =', x_max_ssim, '| y_max_ssim   =',
                  y_max_ssim,   "| r=", np.max(ccrr_ssim))
            return  x_max_ssim/scale_factor, y_max_ssim/scale_factor, posicion_n_0 , posicion_n_f 
        #-----------------------------------------------------------------------------------------------------------
            if method== 'pearson-ssim':
                print('|----Estimating Pearson coeffients ans Estimating Structural Similarity Index  (SSIM) ---- | ')
                print('WARNING: Estimating the similarity between the images using this method is computationally expensive.')
                for i in tqdm(range(0,iris_interp.shape[0]-255,1), position=0, desc="i", leave=False, colour='green', ncols=100):
                    for j in range(0,iris_interp.shape[1]-255,1):
                        matrix=iris_interp[0+i:253+i,0+j:253+j].copy()
                        corr, valor_p= pearsonr(matrix.flatten(), flatten_alma)
                        ssim_none = ssim(matrix, alma_interp,
                                         data_range=iris_interp.max() - iris_interp.min())

                        x=i+alma_interp.shape[1]//2
                        y=j+alma_interp.shape[0]//2

                        ccrr_pearson[x][y]=corr
                        ccrr_ssim[x][y] = ssim_none
                
                x_max_pearson = np.where(ccrr_pearson==np.max(ccrr_pearson))[1][0]
                y_max_pearson = np.where(ccrr_pearson==np.max(ccrr_pearson))[0][0]

                x_max_ssim = np.where(ccrr_ssim == np.max(ccrr_ssim))[1][0]
                y_max_ssim = np.where(ccrr_ssim == np.max(ccrr_ssim))[0][0]
                
                print('________________________________________________________________________')
                print('x_max_pearson=',x_max_pearson  , '| y_max_pearson=', y_max_pearson,"| r=",np.max(ccrr_pearson))
                
                print('________________________________________________________________________')
                print('x_max_ssim   =', x_max_ssim, '| y_max_ssim   =',
                      y_max_ssim,   "| r=", np.max(ccrr_ssim))
                
                x_mean=  (x_max_pearson+x_max_ssim)/(2*scale_factor)
                x_mean=  round(x_mean)
                
                y_mean=  (y_max_pearson+y_max_ssim)/(2*scale_factor)
                y_mean=  round(y_mean)
                return   x_mean/scale_factor, y_mean/scale_factor, int(posicion_n_0), int(posicion_n_f)
        
    



def aligning_IRIS_with_SDO(sdo_file, iris_data, time_vec_iris,time_vec_sdo,time_vec_ALMA ,path_img_crop):
    """
    This function allows us to find a pixel reference in IRIS images to one snapshot(Photogram)
    to do an alignment of IRIS'S image with SDO'S Images
    For this, the SIMM and Pearson method was used to analyze correlations
    --------------------Input----------------------------------
    Parameters:
               sdo_file:      vector with the frames taken by the SDO/AIA. Its format is .npy 
                              and it was created using the program "funtion_to_create_data_SDO_cube.ipynb"     
               iris_data:     FILE .FITS obtained from IRIS database
               
               time_vec_iris: vector where its components are timeutc of every snapshot (photogram) 
                              of IRIS data cube. The Fortmat of time_vec_iris's components is datetime.datetime, for example
                              one of  time_vec_iris's components look like should be 
                              datetime.datetime(2018, 4, 12, 15, 52, 31, 992187) 
                              THIS VECTOR WAS  BUILT USING THE FUNCTION "creating_time_vector"
                              which goes into more detail
                              
               time_vec_sdo:  vector where its components are timeutc of every snapshot (photogram) 
                              of SDO data cube. The Fortmat of time_vec_sdo's components is datetime.datetime, for example
                              one of  time_vec_sdo's components look like should be 
                              datetime.datetime(2018, 4, 12, 15, 52, 31, 992187)
                              THIS VECTOR WAS BUILT WITH THE PROGRAM "funtion_to_create_data_SDO_cube.ipynb"

               time_vec_ALMA: vector where its components are timeutc of every snapshot (photogram) 
                              of SDO data cube. The Fortmat of time_vec_ALMA's components is datetime.datetime, for example
                              one of  time_vec_ALMA's components look like should be 
                              datetime.datetime(2018, 4, 12, 15, 52, 31, 992187)
                              THIS VECTOR TAKEN FROM  SALAT DATA BASE 
                              
               path_img_crop:
    --------------------Input----------------------------------
    Estimated heliprojective coordinates for the center of the 
    IRIS matrix estimated using the SIMM and PEARSON method
    
    """
    url                 = sorted(glob.glob(path_img_crop))
    sdo_resolution      = fits.open(url[0])[0].header['CDELT2']
    iris_resolution     = iris_data[0].header['CDELT1']
    scale_factor        = iris_resolution/ sdo_resolution
    iris_cube           = iris_data[0].data
    print('scale_factor:', scale_factor)
    print('Rescaling IRIS  images')
    
    rescaled_iris       = []
    for i in tqdm(range(len(iris_cube)), colour='green'):
        rescaled_iris.append(rescale(iris_cube[i], scale_factor, order=0,
                                    anti_aliasing=False, preserve_range=True))
    rescaled_iris       = np.array(rescaled_iris)

    print('________________________________________________________________________')
    print('initial shape:'   , sdo_file.shape,
          '| rescaled shape:', rescaled_iris.shape)

    # -------------------Selecting the same windown of ALMA----------------
    posicion_n_0_sdo, posicion_n_f_sdo   = find_nearest(time_vec_sdo , time_vec_ALMA[0], time_vec_ALMA[-1],'SDO', 'ALMA')
    posicion_n_0_iris, posicion_n_f_iris = find_nearest(time_vec_iris, time_vec_ALMA[0], time_vec_ALMA[-1],'IRIS', 'ALMA')
    
    mean_sdo    = np.sum([i for i in sdo_file[int(posicion_n_0_sdo):int(posicion_n_f_sdo)]], 0)/len(sdo_file[int(posicion_n_0_sdo):int(posicion_n_f_sdo)])
    mean_iris   = np.sum([i for i in rescaled_iris[int(posicion_n_0_iris):int(posicion_n_f_iris)]],0)/len(rescaled_iris[int(posicion_n_0_iris):int(posicion_n_f_iris)])

    sdo_interp  = np.interp(mean_sdo, (np.percentile(mean_sdo, 11), np.percentile(mean_sdo, 95)), (0, +1))
    iris_interp = np.interp(mean_iris, (np.percentile(mean_iris, 11), np.percentile(mean_iris, 95)), (0, +1))
    
    
    # -----------------------ALMA_VECTOR---------------------------------------
    flatten_iris = iris_interp.flatten()
    ccrr_pearson = np.zeros((sdo_interp.shape[0], sdo_interp.shape[1]))
    ccrr_ssim    = np.zeros((sdo_interp.shape[0], sdo_interp.shape[1]))
    print('|----Estimating Pearson and SSIM coeffients---- | ')
    for i in tqdm(range(0, sdo_interp.shape[0]-iris_interp.shape[0]+1, 1), position=0, desc="i", leave=False, colour='green', ncols=100):
        for j in range(0, sdo_interp.shape[1] - iris_interp.shape[1]+1, 1):
            matrix = sdo_interp[0+i:iris_interp.shape[0]+i, 0+j:iris_interp.shape[1]+j].copy()
            corr, valor_p = pearsonr(matrix.flatten(), flatten_iris)
            ssim_none = ssim(matrix, iris_interp,
                             data_range=sdo_interp.max() - sdo_interp.min())
            x = i+iris_interp.shape[1]//2
            y = j+iris_interp.shape[0]//2
            ccrr_pearson[x][y] = corr
            ccrr_ssim[x][y]    = ssim_none

    x_max_pearson = np.where(ccrr_pearson == np.max(ccrr_pearson))[1][0]
    y_max_pearson = np.where(ccrr_pearson == np.max(ccrr_pearson))[0][0]
    x_max_ssim    = np.where(ccrr_ssim    == np.max(ccrr_ssim))[1][0]
    y_max_ssim    = np.where(ccrr_ssim    == np.max(ccrr_ssim))[0][0]

    print('________________________________________________________________________')
    print('x_max_pearson=', x_max_pearson, '| y_max_pearson=',
          y_max_pearson, "| r=", np.max(ccrr_pearson))

    
    print('________________________________________________________________________')
    print('x_max_ssim   =', x_max_ssim, '| y_max_ssim   =',
          y_max_ssim,   "| r=", np.max(ccrr_ssim))

    map_sunpy                            = sunpy.map.Map(url[posicion_n_0_sdo])
    x_scaled_pearson , y_scaled_pearson  = x_max_pearson,  y_max_pearson
    cord_pearson                         = map_sunpy.pixel_to_world(x_scaled_pearson *u.pix,  y_scaled_pearson *u.pix )
    print('Pearson: (Tx, Ty) in arcsec',cord_pearson.Tx.value,cord_pearson.Ty.value)
 
    print('________________________________________________________________________')
    x_scaled_ssim, y_scaled_ssim =  x_max_ssim   ,  y_max_ssim
    cord_ssim                    = map_sunpy.pixel_to_world(x_scaled_ssim*u.pix,  y_scaled_ssim*u.pix )
    print('SSIM   : (Tx, Ty) in arcsec',cord_ssim.Tx.value,cord_ssim.Ty.value)






def aligning_IRIS_with_SDO_two(sdo_file, iris_data, iris_resolution,time_vec_iris,time_vec_sdo,time_vec_ALMA ,path_img_crop):
    """
    This function allows us to find a pixel reference in IRIS images to one snapshot(Photogram)
    to do an alignment of IRIS'S image with SDO'S Images
    For this, the SIMM and Pearson method was used to analyze correlations
    --------------------Input----------------------------------
    Parameters:
               sdo_file:      vector with the frames taken by the SDO/AIA. Its format is .npy 
                              and it was created using the program "funtion_to_create_data_SDO_cube.ipynb"     
               iris_data:      data cube obtained from IRIS database in format NP array 
               
               time_vec_iris: vector where its components are timeutc of every snapshot (photogram) 
                              of IRIS data cube. The Fortmat of time_vec_iris's components is datetime.datetime, for example
                              one of  time_vec_iris's components look like should be 
                              datetime.datetime(2018, 4, 12, 15, 52, 31, 992187) 
                              THIS VECTOR WAS  BUILT USING THE FUNCTION "creating_time_vector"
                              which goes into more detail
                              
               time_vec_sdo:  vector where its components are timeutc of every snapshot (photogram) 
                              of SDO data cube. The Fortmat of time_vec_sdo's components is datetime.datetime, for example
                              one of  time_vec_sdo's components look like should be 
                              datetime.datetime(2018, 4, 12, 15, 52, 31, 992187)
                              THIS VECTOR WAS BUILT WITH THE PROGRAM "funtion_to_create_data_SDO_cube.ipynb"

               time_vec_ALMA: vector where its components are timeutc of every snapshot (photogram) 
                              of SDO data cube. The Fortmat of time_vec_ALMA's components is datetime.datetime, for example
                              one of  time_vec_ALMA's components look like should be 
                              datetime.datetime(2018, 4, 12, 15, 52, 31, 992187)
                              THIS VECTOR TAKEN FROM  SALAT DATA BASE 
                              
               path_img_crop:
    --------------------Input----------------------------------
    Estimated heliprojective coordinates for the center of the 
    IRIS matrix estimated using the SIMM and PEARSON method
    
    """


    url                 = sorted(glob.glob(path_img_crop))
    sdo_resolution      = fits.open(url[0])[0].header['CDELT2']
    scale_factor        = iris_resolution/ sdo_resolution
    iris_cube           = iris_data
    print('scale_factor:', scale_factor)
    print('Rescaling IRIS  images')
    
    rescaled_iris       = []
    for i in tqdm(range(len(iris_cube)), colour='green'):
        rescaled_iris.append(rescale(iris_cube[i], scale_factor, order=0,
                                    anti_aliasing=False, preserve_range=True))
    rescaled_iris       = np.array(rescaled_iris)

    print('________________________________________________________________________')
    print('initial shape:'   , sdo_file.shape,
          '| rescaled shape:', rescaled_iris.shape)

    # -------------------Selecting the same windown of ALMA----------------
    posicion_n_0_sdo, posicion_n_f_sdo   = find_nearest(time_vec_sdo , time_vec_ALMA[0], time_vec_ALMA[-1],'SDO', 'ALMA')
    posicion_n_0_iris, posicion_n_f_iris = find_nearest(time_vec_iris, time_vec_ALMA[0], time_vec_ALMA[-1],'IRIS', 'ALMA')
    
    mean_sdo    = np.sum([i for i in sdo_file[int(posicion_n_0_sdo):int(posicion_n_f_sdo)]], 0)/len(sdo_file[int(posicion_n_0_sdo):int(posicion_n_f_sdo)])
    mean_iris   = np.sum([i for i in rescaled_iris[int(posicion_n_0_iris):int(posicion_n_f_iris)]],0)/len(rescaled_iris[int(posicion_n_0_iris):int(posicion_n_f_iris)])

    sdo_interp  = np.interp(mean_sdo, (np.percentile(mean_sdo, 11), np.percentile(mean_sdo, 95)), (0, +1))
    iris_interp = np.interp(mean_iris, (np.percentile(mean_iris, 11), np.percentile(mean_iris, 95)), (0, +1))
    
    
    # -----------------------ALMA_VECTOR---------------------------------------
    flatten_iris = iris_interp.flatten()
    ccrr_pearson = np.zeros((sdo_interp.shape[0], sdo_interp.shape[1]))
    ccrr_ssim    = np.zeros((sdo_interp.shape[0], sdo_interp.shape[1]))
    print('|----Estimating Pearson and SSIM coeffients---- | ')
    for i in tqdm(range(0, sdo_interp.shape[0]-iris_interp.shape[0]+1, 1), position=0, desc="i", leave=False, colour='green', ncols=100):
        for j in range(0, sdo_interp.shape[1] - iris_interp.shape[1]+1, 1):
            matrix = sdo_interp[0+i:iris_interp.shape[0]+i, 0+j:iris_interp.shape[1]+j].copy()
            corr, valor_p = pearsonr(matrix.flatten(), flatten_iris)
            ssim_none = ssim(matrix, iris_interp,
                             data_range=sdo_interp.max() - sdo_interp.min())
            x = i+iris_interp.shape[1]//2
            y = j+iris_interp.shape[0]//2
            ccrr_pearson[x][y] = corr
            ccrr_ssim[x][y]    = ssim_none

    x_max_pearson = np.where(ccrr_pearson == np.max(ccrr_pearson))[1][0]
    y_max_pearson = np.where(ccrr_pearson == np.max(ccrr_pearson))[0][0]
    x_max_ssim    = np.where(ccrr_ssim    == np.max(ccrr_ssim))[1][0]
    y_max_ssim    = np.where(ccrr_ssim    == np.max(ccrr_ssim))[0][0]

    print('________________________________________________________________________')
    print('x_max_pearson=', x_max_pearson, '| y_max_pearson=',
          y_max_pearson, "| r=", np.max(ccrr_pearson))

    
    print('________________________________________________________________________')
    print('x_max_ssim   =', x_max_ssim, '| y_max_ssim   =',
          y_max_ssim,   "| r=", np.max(ccrr_ssim))

    map_sunpy                            = sunpy.map.Map(url[posicion_n_0_sdo])
    #x_scaled_pearson , y_scaled_pearson  = x_max_pearson,  y_max_pearson
    cord_pearson                         = map_sunpy.pixel_to_world(x_max_pearson *u.pix,  y_max_pearson *u.pix )
    print('Pearson: (Tx, Ty) in arcsec',cord_pearson.Tx.value,cord_pearson.Ty.value)
 
    print('________________________________________________________________________')
    #x_scaled_ssim, y_scaled_ssim =  x_max_ssim   ,  y_max_ssim
    cord_ssim                    = map_sunpy.pixel_to_world(x_max_ssim*u.pix,  y_max_ssim*u.pix )
    print('SSIM   : (Tx, Ty) in arcsec',cord_ssim.Tx.value,cord_ssim.Ty.value)


        
    x_cord_hep = (cord_ssim.Tx.value+cord_pearson.Tx.value)/2
    y_cord_hep = (cord_ssim.Ty.value+cord_pearson.Ty.value)/2
    
    return x_cord_hep, y_cord_hep
    
 
def sdo_rescale(sdo_file, sdo_resolution,alma_resolution):
    '''
    This function allows us to rescale the SDO images so that the 
    pixel size of the SDO is the same as that of ALMA
    --------------------Input----------------------------------
    Parameters:  
              sdo_file: type array   is a vector where its position correspond a matrices 
                           (A photogram of  the sum)  taken from SDO 
    
    Output:   SDO Images rescaled (a datacube)
    '''
    scale_factor = sdo_resolution/alma_resolution
    print('scale_factor:', scale_factor)
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
    return rescaled_sdo



def cutting_sdo_data(cube_data,cube_time,files,cord_aling, alma_cube,alma_time ,wave_len ,alma_resolution):
    aperture           = aperture_circ_ALMA(alma_cube)
    sdo_resolution     = fits.open(files[0])[0].header['CDELT1']
    rescaling_cubo     = sdo_rescale(cube_data ,sdo_resolution , alma_resolution)
    n_0_time, n_f_time = find_nearest(cube_time, alma_time[0], alma_time[-1], 'AIA {}'.format(wave_len), 'ALMA')
    files=files[n_0_time]
    aia_rot     = sunpy.map.Map(files)
    scale_sdo   = sdo_resolution/alma_resolution
    c_heliopro      = SkyCoord(Tx=cord_aling[0]* u.arcsec  , Ty=cord_aling[1]*u.arcsec, frame=aia_rot.coordinate_frame)
    x_max_sdo, y_max_sdo= aia_rot.world_to_pixel(c_heliopro).x.value*scale_sdo, aia_rot.world_to_pixel(c_heliopro).y.value*scale_sdo
    
    
    print(x_max_sdo, y_max_sdo,scale_sdo )
    y_new_i_sdo=int(y_max_sdo-alma_cube.shape[1]/2)
    y_new_f_sdo=int(y_max_sdo+alma_cube.shape[1]/2)
    x_new_i_sdo=int(x_max_sdo-alma_cube.shape[2]/2)
    x_new_f_sdo=int(x_max_sdo+alma_cube.shape[2]/2)

    if  int(y_new_f_sdo- y_new_i_sdo)!=alma_cube.shape[1]:
        y_new_f_sdo=alma_cube.shape[1]+y_new_i_sdo

    if  int(x_new_f_sdo- x_new_i_sdo)!=alma_cube.shape[2]:
        x_new_f_sdo=alma_cube.shape[2]+x_new_i_sdo   

    image_cuted=[]
    for i in range(n_0_time,n_f_time+1,1):
        image_cuted.append(rescaling_cubo[i][y_new_i_sdo:y_new_f_sdo, x_new_i_sdo:x_new_f_sdo]*aperture)
    image_cuted=np.array(image_cuted)
    cube_time=cube_time[n_0_time:n_f_time+1]
    
    plt.figure()              
    plt.title('WAVELNTH {} Å'.format(fits.open(files)[0].header['WAVELNTH']))
    vmin=np.percentile(image_cuted,11)
    vmax=np.percentile(image_cuted,95)
    plt.imshow(image_cuted[0],  origin='lower', vmin=vmin , vmax=vmax)
    plt.colorbar()
    plt.show()
    return image_cuted, cube_time
    
    
   
    
def  cutting_iris_data(files, time_vec,cord_pix_max  ,alma_cube, alma_time , pix_size_alma):
    aperture           = aperture_circ_ALMA(alma_cube)
    rescaling_iris     = iris_rescale(files, pix_size_alma)
    aperture           = aperture_circ_ALMA(alma_cube)
    scale_factor       =  files[0].header['CDELT1']/pix_size_alma
    n_0_time, n_f_time = find_nearest(time_vec, alma_time[0], alma_time[-1], 'IRIS', 'ALMA')
    time_vec_iris = creating_time_vector(files)
    x_resul_iris=cord_pix_max[0]*scale_factor
    y_resul_iris=cord_pix_max[1]*scale_factor
    radio_new_pix= alma_cube.shape[1]/2

    x_0=int(x_resul_iris-radio_new_pix)
    x_f=int(x_resul_iris+radio_new_pix)
    y_0=int(y_resul_iris-radio_new_pix)
    y_f=int(y_resul_iris+radio_new_pix)
    
    if  int(x_f- x_0)!=alma_cube.shape[2]:
        x_f=alma_cube.shape[2]+x_0
    if  int(y_f- y_0)!=alma_cube.shape[1]:
        y_f=alma_cube.shape[1]+y_0

    image_cuted=[]
    for i in range(n_0_time,n_f_time+1,1):
        image_cuted.append(rescaling_iris[i][y_0:y_f, x_0:x_f]*aperture)
    image_cuted=np.array(image_cuted)
    time_vec=time_vec[n_0_time:n_f_time+1]
    
    plt.figure()
    plt.title(str(time_vec_iris[n_0_time]))
    vmin= np.percentile(image_cuted,13) 
    vmax= np.percentile(image_cuted,98) 
    plt.imshow(image_cuted[0],  origin='lower', vmin=vmin , vmax=vmax)
    plt.colorbar()
    plt.show()
    
    np.save('result_align/crop_and_aling_iris_{}_Pixel_ALMA_size'.format(int(files[0].header['TWAVE1'])), image_cuted , allow_pickle=True)
    np.save('result_align/time_cut_iris_{}'.format(int(files[0].header['TWAVE1'])), time_vec , allow_pickle=True)
    
    
    
    
    
    
def find_nearest_modificated(array, alma_t0):
    '''
    This function allows us to find (in a time vector) 
    the initial and final times at times alma_t0, alma_tf.
    
    --------------------Input----------------------------------
    Parameters: 
                array : type array and its elements is in the datetime.datetime 
                        format (in this case we  obtained this vector using the )
                        in our case this time vector was obtained with the DATA-OBS 
                        of each previously downloaded SDO frame  
                alma_t0: type T0 is  where the ALMA observation has started
                         datetime.datetime  format
                alma_tf: type Tf is  where the ALMA observation has finalized 
                         datetime.datetime  format  
    --------------------Output----------------------------------
    Positions in the vector array (array) closest to the   alma_t0, alma_tf        
    '''
    array = np.asarray(array)
    idx_0 = (np.abs(array - alma_t0)).argmin()
    #print(array[idx_0],alma_t0)
    return idx_0
