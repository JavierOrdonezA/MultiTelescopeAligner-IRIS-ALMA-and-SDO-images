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





![Alineación ALMA con IRIS](IMAGES_RESULTS/result_alingnment_with_iris.jpg)




![Alineación ALMA con IRIS](IMAGES_RESULTS/IRIS_alignment_with_SDO.png)



