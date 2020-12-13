# Region-Filling-and-Object-Removal-by-Exemplar-Based-Image-Inpainting

This repository contains the implementation of [Region-Filling-and-Object-Removal-by-Exemplar-Based-Image-Inpainting](https://www.irisa.fr/vista/Papers/2004_ip_criminisi.pdf) in MATLAB. This was implemented as a part of assingment in Digital Image Processing Course, E & ECE, IITKGP.

Team Members:
1. Sirish Gambhira <br>
2. Sanku Yogesh <br>
3. MSR Pavan Kumar <br>
4. Bandaru Sai Manoj <br>

## Requirements
1. Curve fitting
2. Basic plot tool boxes

## Code outline

For the difference between Code1, Code2 refer 6th point of discussion <br>

AssistedFreehand is used to get the mask, hold the left button to draw, and double click when completed <br>

Functions:<br>
1. In the first two blocks of code, we are reading the image and taking the user interactive mask.
2. Get-contour computes the contour where the bright pixels in the binary image with any of the 8-neighboring pixels as dark are considered in the contour
3. get_normal computes the normal on each point of contour by using the curve fitting method where we arrange the points in the clockwise direction and interpolated the intermediate values.
4. get_priorities computes the confidence and data values for each pixel on the contour and techniques used are discussed in the report.
5. get_max computes the index of the pixel having the maximum priority values.
6. get_sample and get_distance functions search the best matching exemplar.
7. update function updates the Image, binary image/mask, and confidence values
8. Check function checks whether the task is completed or not

## Results





