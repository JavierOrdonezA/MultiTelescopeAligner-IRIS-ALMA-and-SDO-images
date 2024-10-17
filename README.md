# Alignment of IRIS, ALMA, and SDO Images
  
## Introduction  
This project develops a  method for aligning images from three major solar observation instruments: the Atacama Large Millimeter/submillimeter Array (ALMA), Interface Region Imaging Spectrograph (IRIS), and the Solar Dynamics Observatory (SDO). The code focuses on co-observations in the D06 region within the Solar ALMA Science Archive (SALSA), enhancing the understanding of solar atmospheric dynamics by synchronizing data from these telescopes spatially and temporally.

## Motivation and Objectives
The dynamic nature of the solar atmosphere presents significant challenges in solar physics research, particularly when comparing data across different observation platforms. The primary objective of this project is to address the lack of a dedicated preprocessing technique for the complex analysis of multispectral solar data, thereby improving the quality and reliability of data comparisons.

## Features
- Precise alignment leveraging SDO helioprojective coordinates.
- Utilizes Pearson Correlation Coefficient (PCC) and Structural Similarity Index (SSIM) for robust image comparison.
- Applicable to diverse observational data from ALMA, IRIS, and SDO.

## Installation
To install this library, please follow these steps:

1. Clone the repository from GitHub:
   ```bash
   git clone https://github.com/JavierOrdonezA/Alignment-of-IRIS-ALMA-and-SDO-images.git
2. Navigate to the cloned directory:
   ```bash
   cd Alignment-of-IRIS-ALMA-and-SDO-images
3. Execute the setup script to install the library:
   ```bash
   python setup.py install
   
This will install all necessary dependencies and make the library available for use in your Python environment.

## Usage and Examples
Detailed examples on how to utilize the functions are provided in the `Guide_Example_and_Results` folder. This includes a comprehensive guide on the types of data (formats) that should be used as input for the functions. Please refer to the `Guide_Example_and_Results` directory to understand the setup and execution of functions within this library.

## Methodology

### Data Acquisition

Data is systematically sourced from ALMA, IRIS, and SDO, focusing on regions observed simultaneously by these telescopes. The SALSA database is instrumental in identifying relevant observational data.

### Data Preparation
Images from SDO and ALMA are rotated to align with solar north. IRIS images are rescaled to match ALMA’s pixel size using the arcsecond ratio.

### Alignment Process
* ALMA and IRIS: Correlation analysis is performed between processed images to find the best alignment using PCC and SSIM.
* IRIS and SDO: Cropped and rescaled IRIS images are aligned with corresponding SDO frames, enhancing the Field of View to slightly exceed that of ALMA’s.

## Results and Analysis
Figures included in `Guide_Example_and_Results` demonstrate the successful correlation and alignment of the thermal structures observed by ALMA with data from IRIS and SDO.

![Alignment Results](https://github.com/JavierOrdonezA/MultiTelescopeAligner-IRIS-ALMA-and-SDO-images/blob/master/Guide_Example_and_Results/result_alingnment_with_iris_new_plot_dpi_100.jpg)


## Conclusions
This alignment methodology not only facilitates a deeper understanding of the solar atmosphere’s dynamics but also proves highly accurate and computationally efficient. The approach is adaptable and can potentially be extended to other telescopes and observational data sets.


## Future Directions
Future updates may include an advanced alignment algorithm that considers projection effects, aiming to further refine the precision and applicability of this technique.

## License
This project is open source and available under the MIT License. This allows for modification, distribution, and private use.

MIT License

Copyright (c) [2024] [F. J. Ordonez Araujo, J. C Guevara Gomez]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Library"), to deal
in the Library without restriction, including without limitation the rights
to use, copy, modify, merge, distribute, and sublicense.




## Citation

For more details on the process, see:
**Method article**: 'https://arxiv.org/abs/2404.04401',

**Second chapter of my master's thesis** for more details than in the paper: 'https://repositorio.unal.edu.co/handle/unal/85838'.

---
Feel free to explore the projects and their solutions in the respective directories.
👾 Happy coding! 🥷

**F. J. Ordonez Araujo (fordonezaraujo@gmail.com)**

**J. C Guevara Gomez (juancamilo.guevaragomez@gmail.com)**

---
Thanks For Watching This Repository!

**KEEP AWESOME & STAY COOL!** 😎

Feel Free To Fork And Report If You Find Any Issue :) 




