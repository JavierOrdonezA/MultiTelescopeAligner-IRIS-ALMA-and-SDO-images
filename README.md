# Alignment of IRIS, ALMA, and SDO Images

## Introduction
This project develops a robust method for aligning images from three major solar observation instruments: the Atacama Large Millimeter/submillimeter Array (ALMA), Interface Region Imaging Spectrograph (IRIS), and the Solar Dynamics Observatory (SDO). The code focuses on co-observations in the D06 region within the Solar ALMA Science Archive (SALSA), enhancing the understanding of solar atmospheric dynamics by synchronizing data from these telescopes spatially and temporally.

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
