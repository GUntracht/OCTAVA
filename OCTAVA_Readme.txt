*************** OCTAVA Help ***************

This file contains help and other data for the 
OCTAVA vascular analysis toolbox 

ImageJ must be launched from within the OCTAVA app to run the vascular analysis

*******************************************
*******************************************
System Requirements

MATLAB APP version: 
MATLAB version must be 2019 or more recent
MATLAB toolboxes required: 
-Fuzzy Logic toolbox 

Standalone version: 
MATLAB runtime environment 2020b 

*******************************************
*******************************************
Installation instructions 
For a pre-compiled standalone version, please contact gavuntracht at Gmail dot com. The pre-compiled version can be installed automatically by running the installation file that we will provide. 
For the GitHub version:
1.	Download OCTAVA.mlprj and the folder scripts into your local directory
2.	Copy the ImageJ/FIJI installation and source files into your directory 
	If you are downloading the latest version of ImageJ, you may need to update the copoytoImgPlus macro. 
	update copytoImgPlus in path .../ImageJ/Fiji/scripts by pulling: https://github.com/fiji/fiji/blob/master/scripts/copytoImgPlus.m
3.	Save mij.jar and imagej-matlab-0.7.2 in the ImageJ folder
4.	Save OCTAVA.ijm into the folder ...ImageJ\Fiji\macros\OCTAVA Tools. You may need to create the folder named ¨OCTAVA Tools¨.

*******************************************
*******************************************

Process single image instructions

1. Load image into OCTAVA
2. Binarize to generate binary mask 
3. Frangi vesselness filter can be turned on and off, or sigma_max can 
	be adjusted to accommodate larger or smaller maximum vessel diameters 
4. Press "skeletonize" to skeletonize, identify vessels and nodes
5. Press "Generate Metrics" to calculate network metrics and generate histograms 
    See below for more details 
   

Resize
For large images, the ImageJ algorithm may run out of memory depending on 
the performance of your computer. Resize gives you the option to compress 
images to reduce the number of pixels. This may impact the precision of your
results 

Image subset
Allows the user to select a rectangular or circular region from within the 
loaded image. Select the desired shape and press the 'Start' button. The 
ROI can then be drawn onto the image directly. 

Median Filter
Removes small noise features from the image size using a smoothing filter
The numerical value represents the filter size in pixels. A larger number 
corresponds to removal of larger features

Binarization Method
The user can select between two binarization algorithms: Fuzzy Means and 
Adaptive Thresholding. 

Frangi Filter
The Frangi filter is a Hessian filter which enhances vessel-like structures 
in the image. Frangi S_Max gives an indication of the upper bound of vessel 
diameter expected in the image. A larger number corresponds to a larger 
maximum vessel diameter. Overestimation of Frangi S_max can lead to 
artificial vessel dilation. 

Invert
Binary mask must be generated with the vessels in white and background in
black. If the binary mask is generated with black vessels, the Invert button
can be use to switch the color values of the two regions. 

Manual Curation
Checking this box allows the user to manually add and remove ROIs from processing
To add branches, draw a line on the image and press 't'
To remove, delete the relevant ROI from the ImageJ ROI manager 
Only master segments, branches, isolated elements, and junctions are included 
in the metric calculations, so these should be the focus when removing artefacts. 
Aby added ROIs are counted as branches. 

Twig Size
This setting determines the maximum diameter of the of twigs (in pixels)
Isolated twigs are excluded from statistical analysis

Generate Metrics
The size of the image in mm must be specified in order to obtain accurate 
results. The size here should refer to the original loaded image before 
cropping or resizing. If the image is circular, the diameter of the image 
should be specified. 



********************************************

Saving data:
From the Saving menu:

-Save Binary Map - saves only an image of the binary map 
-Save Metrics as mat File - saves the metrics and histogram data from the currently active 
    image in a .mat file 
-Save Histogram Data - saves the list of diameters, tortuosity, and lengths 
    from the currently active image in an Excel file  
-Save to Excel - saves a record of metrics  for all datasets processed within the session 
-Clear Session - clears the record of datasets processed within the session 

Images of the histograms can also be saved directly by hovering the mouse 
over the title of the plot. The icon with the arrow gives options for export 
including saving (disc icon) or copying directly. 

ImageJ menu:
If needed, this menu can be used to open and close ImageJ from within MOVA

OCTAVA menu:
Clear Plots - clears and resets all plots in the app

********************************************

Batch mode instructions 
1. Select the location of raw OCT images and specify saving directory 
2. Press start batch process to process all images in the specified folder 
3. Biomarkers will be automatically generated and saved in an Excel spreadsheet 

The file type must be selected for OCTAVA to locate the images withing the
selected folder. 

The user can pre-binarize images in order to optimize the binarization 
individually for different images and still use batch process for generation 
of metrics. In this case, the overlay images will be saved on the binary
 mask 
instead of the OCT image since OCTAVA will not have access to the original OCT images. 


Use esc to abort ImageJ Macro



********************************************
Visualizing plots and images 

Hovering the mouse over the title of the plot opens a plot manipulation 
window.

Save icon: allows the user to directly save the plots from the GUI
Magnifying glasses: allow the user to zoom in and out 
Home icon: returns the plot to the original scaling and position centered
around the image or data 

********************************************
Rights and Permissions

OCTAVA employs the following existing software packages

1. ImageJ (publicly available) 
 ImageJ (http://rsb.info.nih.gov/ij/) is a public domain Java image processing program inspired by NIH Image for the Macintosh.
It runs, either as an online applet or as a downloadable application, on any computer with a Java 1.1 or later virtual machine. 
Downloadable distributions are available for Windows, Mac OS, Mac OS X and Linux. 
The author, Wayne Rasband (wayne@codon.nih.gov), is at the Research Services Branch, National Institute of Mental Health, 
Bethesda, Maryland, USA.
 
2. Angiogenesis analyzer in ImageJ  
Written by Gilles Carpentier, 2008. The macro is available at http://rsb.info.nih.gov/ij/macros/toolsets/Dot%20Blot%20Analyzer.txt 
and more information can be found at http://image.bio.methods.free.fr/dotblot.html). 

3. ImageJ-MATLAB package 
Freely available extension to ImageJ2 (http://imagej.net/Downloads). Installation and use instructions available at 
http://imagej.net/MATLAB_Scripting. Tested with ImageJ 2.0.0-rc-54, Java 1.8.0_66 and MATLAB R2015b. 
More information is available at https://doi.org/10.1093/bioinformatics/btw681.
