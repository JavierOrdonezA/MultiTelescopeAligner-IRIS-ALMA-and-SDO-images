# Alignment-of-IRIS-ALMA-and-SDO-images
This project presents a code for aligning images from ALMA, IRIS, and SDO telescopes with co-observations. It includes functions for aligning IRIS with ALMA, IRIS with SDO, and SDO with ALMA. The code aligns images from the D06 region of the Solar ALMA Science Archive (SALSA) database.

## Authors and Acknowledgments
- Francisco J. Ordoñez Araujo
- Juan Camilo Guevara Gómez
- Benjamín Calvo Mozo


## Motivation and objetives
The primary motivation for this project arises from the challenges associated with studying the dynamic and continuously changing physical conditions of the solar atmosphere. Accurately comparing different regions and understanding the underlying physical phenomena require robust data processing techniques. 

In the field of solar physics, research often emphasizes specific research objectives, leading to limited discussions about data preprocessing in academic literature. This gap highlights the need for a dedicated preprocessing method that can address the complexities of analyzing data from telescopes like IRIS, SDO, and ALMA.

Recognizing this need, our project introduces an innovative preprocessing method designed to spatially and temporally synchronize images from these diverse telescopes. This alignment technique is crucial for enhancing the quality and reliability of data analysis, thereby enabling a deeper and more accurate understanding of the solar atmosphere's dynamics.


## Features
- Alignment based on helioprojective coordinates.
- Use of Pearson Correlation Coefficient (PCC) and Structural Similarity Index (SSIM).
- Application to solar observations from ALMA, IRIS, and SDO.

## Method Used
The methodology involves a series of steps designed to co-align solar observations from different telescopes. Using the Solar Dynamics Observatory (SDO) helioprojective coordinates as a reference, we align ALMA and IRIS images. This process incorporates both Pearson Correlation Coefficient and Structural Similarity Index analysis to establish the highest degree of correlation between the datasets.






## Methodology

### Overview
This project involves a comprehensive approach to align solar images captured by various telescopes: Helioseismic and Magnetic Image (HMI), Atmospheric Imaging Assembly (AIA) on the Solar Dynamics Observatory (SDO), Interface Region Imaging Spectrograph (IRIS), and Atacama Large Millimeter/submillimeter Array (ALMA). The focus is on region D06, observed on April 12, 2018, with data encompassing temperature, brightness, and magnetic field measurements.

### Alignment Process
The alignment process is essential for correcting deviations in helioprojective coordinates due to factors like atmospheric turbulence and instrumental errors. It comprises four main stages:

1. **Data Acquisition:** 
   - Data from ALMA, IRIS, and SDO is precisely downloaded and selected, focusing on observations of the same solar region.
   - The SALSA database is utilized for identifying regions observed by ALMA, which are co-observed with IRIS and SDO.

2. **Data Preparation:**
   - Images from SDO and ALMA are rotated for proper alignment, using specific header values to match solar north.

3. **Alignment of ALMA Images with IRIS:**
   - IRIS images are rescaled to match ALMA's pixel size.
   - Averaging and interpolating IRIS and ALMA images to create two comparable datasets.
   - Correlation analysis using Pearson Correlation Coefficient (PCC) and Structural Similarity Index (SSIM) to align images.

4. **Alignment of IRIS Images with SDO:**
   - Cropping and normalizing IRIS images, followed by rescaling to match SDO's pixel size.
   - Correlation analysis to align IRIS images with SDO, thereby aligning ALMA images with SDO by transitivity.

### Optimization and Implementation
The process is optimized to limit computational costs, including strategic cropping to narrow the search range for maximum correlations. The entire procedure is implemented in Python, More details of this process will be found in an article that will be published, whose DOI will be added to this README when it is published.





# Alignment of IRIS, ALMA, and SDO Images

## Introduction
This project presents a  code for aligning images from the Atacama Large Millimeter/submillimeter Array (ALMA), Interface Region Imaging Spectrograph(IRIS), and Solar Dynamics Observatory (SDO) telescopes, focusing on co-observations. The code specifically targets the alignment of images from the D06 region within the Solar ALMA Science Archive (SALSA) database, offering functions to align IRIS with  ALMA, IRIS with SDO, and SDO with ALMA.

## Motivation and objetives
The primary motivation for this project arises from the challenges associated with studying the dynamic and continuously changing physical conditions of the solar atmosphere. Accurately comparing different regions and understanding the underlying physical phenomena require robust data processing techniques. 

In the field of solar physics, research often emphasizes specific research objectives, leading to limited discussions about data preprocessing in academic literature. This gap highlights the need for a dedicated preprocessing method that can address the complexities of analyzing data from telescopes like IRIS, SDO, and ALMA.

Recognizing this need, our project introduces an innovative preprocessing method designed to spatially and temporally synchronize images from these diverse telescopes. This alignment technique is crucial for enhancing the quality and reliability of data analysis, thereby enabling a deeper and more accurate understanding of the solar atmosphere's dynamics.

## Features
- Precise alignment based on SDO helioprojective coordinates.
- Utilization of Pearson Correlation Coefficient (PCC) and Structural Similarity Index (SSIM) for image comparison.
- Application to diverse solar observations from ALMA, IRIS, and SDO.

## Methodology

### Data Acquisition
- Systematic download and selection of data from ALMA, IRIS, and SDO, with a focus on observations of the same solar region.
- Utilization of the SALSA database for identifying regions observed by ALMA in co-observation with IRIS and SDO.

### Data Preparation
- Rotation of SDO and ALMA images using header values to align with solar north.

### Alignment Process
1. **ALMA Images with IRIS:**
   - Rescaling of IRIS images to match ALMA's pixel size.
   - Averaging and interpolating images for comparability.
   - Correlation analysis using PCC and SSIM.

2. **IRIS Images with SDO:**
   - Cropping and normalization of IRIS images, followed by rescaling to match SDO's pixel size.
   - Correlation analysis for alignment.

### Optimization and Implementation
The process is optimized to limit computational costs, including strategic cropping to narrow the search range for maximum correlations. The entire procedure is implemented in Python.

### Upcoming Publication
More details on this process will be available in an upcoming article, with the DOI to be added to this README upon publication.


## Authors and Acknowledgments
- Francisco J. Ordoñez Araujo
- Juan Camilo Guevara Gómez
- Benjamín Calvo Mozo

We extend our gratitude to all collaborators and contributors to this project.









![Alineación ALMA con IRIS](IMAGES_RESULTS/result_alingnment_with_iris.jpg)




![Alineación ALMA con IRIS](IMAGES_RESULTS/IRIS_alignment_with_SDO.png)



