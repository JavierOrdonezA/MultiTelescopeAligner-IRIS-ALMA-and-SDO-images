import unittest
from unittest.mock import MagicMock, patch

from astropy.coordinates import SkyCoord
from skimage.transform import rescale
from datetime import  date, datetime,timedelta
import matplotlib
from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np 
from skimage.metrics import structural_similarity as compare_ssim
from sunpy.net import Fido, attrs as SAF
from  scipy.stats.stats import pearsonr
from skimage.metrics import structural_similarity as ssim
import scipy.ndimage as scnd
import sunpy.map
import numpy as np 
import sys
import sunpy
import glob 
import os
import warnings



sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../solar_alignment_py')))

from myfunctions import (creating_time_vector,
			fun_alma_rotate, 
			aperture_circ_ALMA,
			find_nearest,
			create_data_cube,
			calculate_correlation,
			calculate_correlation_pearson,
			calculate_ssim,
			get_query)
			
warnings.filterwarnings("ignore")



class TestCreatingTimeVector(unittest.TestCase):

    def setUp(self):
        # Mocking the iris_file input
        self.mock_iris_file = MagicMock()
        self.mock_iris_file[0].header = {
            'CADPL_AV': 12,  # Cadence in seconds
            'STARTOBS': '2023-01-01T00:00:00.000'  # Start time
        }
        # Mocking the data array length
        self.mock_iris_file[0].data = np.zeros((100, 100))  # 100 frames

    def test_creating_time_vector(self):
        # Expected outcomes
        expected_start_time = datetime.strptime(self.mock_iris_file[0].header['STARTOBS'], "%Y-%m-%dT%H:%M:%S.%f")
        expected_cadence = self.mock_iris_file[0].header['CADPL_AV']
        
        # Running the function
        time_vector = creating_time_vector(self.mock_iris_file)
        
        # Check type of the output
        self.assertIsInstance(time_vector, np.ndarray, "The output should be a numpy array")
        
        # Check the length of the time_vector
        self.assertEqual(len(time_vector), 100, "The time vector should have 100 elements")
        
        # Check the content and correctness of the time_vector
        for i, time_point in enumerate(time_vector):
            expected_time = expected_start_time + timedelta(seconds=i * expected_cadence)
            self.assertEqual(time_point, expected_time, f"Time vector at index {i} is incorrect")


class TestApertureCircAlma(unittest.TestCase):

    def test_aperture_circ_ALMA(self):
        alma_cube = np.zeros((100, 100, 100))  # Simulating an ALMA cube
        result = aperture_circ_ALMA(alma_cube)
        self.assertEqual(result.shape, (100, 100), "The mask should match the shape of an image in the ALMA cube")
        # Check if the center and a corner point are correct
        self.assertEqual(result[50, 50], 1.0, "Center should be within the aperture")
        self.assertEqual(result[0, 0], 0.0, "Corner should be outside the aperture")


class TestFindNearest(unittest.TestCase):

    def setUp(self):
        # Creating a time array for the test
        start_time = datetime.now()
        self.time_array = [start_time + timedelta(minutes=i) for i in range(60)]  # 60 minutes of data

    def test_find_nearest(self):
        alma_t0 = datetime.now() + timedelta(minutes=20)
        alma_tf = datetime.now() + timedelta(minutes=40)
        idx_0, idx_f = find_nearest(self.time_array, alma_t0, alma_tf, "SDO", "ALMA")
        # Check the indices and the correctness of the found times
        self.assertTrue(abs((self.time_array[idx_0] - alma_t0).total_seconds()) < 60, "Index for t0 should be the closest time")
        self.assertTrue(abs((self.time_array[idx_f] - alma_tf).total_seconds()) < 60, "Index for tf should be the closest time")



class TestCalculateCorrelationPearson(unittest.TestCase):

    def setUp(self):
        # Create sample SDO and IRIS images
        self.sdo_img = np.random.rand(200, 200)  # Larger SDO image
        self.iris_img = np.random.rand(50, 50)  # Smaller IRIS image

    @patch('tqdm.tqdm')  # Mocking tqdm to not affect our tests
    def test_calculate_correlation_pearson(self, mock_tqdm):
        correlation_matrix = calculate_correlation_pearson(self.sdo_img, self.iris_img)

        # Check that the output is a numpy array
        self.assertIsInstance(correlation_matrix, np.ndarray, "Output should be a numpy array")

        # Check the shape of the correlation matrix
        expected_shape = (self.sdo_img.shape[0], self.sdo_img.shape[1])
        self.assertEqual(correlation_matrix.shape, expected_shape, "Correlation matrix should match SDO image dimensions")


class TestCalculateSSIM(unittest.TestCase):

    def setUp(self):
        # Create sample SDO and IRIS images
        self.sdo_img = np.random.rand(200, 200)  # Larger SDO image
        self.iris_img = np.random.rand(50, 50)  # Smaller IRIS image

    @patch('tqdm.tqdm')  # Mocking tqdm to not affect our tests
    def test_calculate_ssim(self, mock_tqdm):
        ssim_matrix = calculate_ssim(self.sdo_img, self.iris_img)

        # Check that the output is a numpy array
        self.assertIsInstance(ssim_matrix, np.ndarray, "Output should be a numpy array")

        # Check the shape of the SSIM matrix
        expected_shape = (self.sdo_img.shape[0], self.sdo_img.shape[1])
        self.assertEqual(ssim_matrix.shape, expected_shape, "SSIM matrix should match SDO image dimensions")




if __name__ == '__main__':
    unittest.main()

