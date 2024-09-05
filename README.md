# OCTAVA

This respoitory contains all the sourcecode required to run the OCTAVA software 
OCTAVA is available free for use under an MIT licence. If you use this software in your research, please cite the following papers: 

Untracht, G.R., Durkee, M.S., Zhao, M. et al. Towards standardising retinal OCT angiography image analysis with open-source toolbox OCTAVA. Sci Rep 14, 5979 (2024). https://doi.org/10.1038/s41598-024-53501-6

Untracht GR, Matos RS, Dikaios N, Bapir M, Durrani AK, Butsabong T, et al. (2021) OCTAVA: An open-source toolbox for quantitative analysis of optical coherence tomography angiography images. PLoS ONE 16(12): e0261052. https://doi.org/10.1371/journal.pone.0261052

This software was last updated in January 2024

In order to install the latest version of the software for editing, please follow the following instructions: 

1. Downoald OCTAVA_2p0.mlapp and the folder scripts into your local directory




For the original version: 

1. Downoald OCTAVA.mlapp and the folder scripts into your local directory
2. Copy the ImageJ/FIJI installation and source files into your directory 
    If you are downloading the latest version of ImageJ, you may need to update the copoytoImgPlus macro.
    update copytoImgPlus
    in path .../ImageJ/Fiji/scripts by pulling:
    https://github.com/fiji/fiji/blob/master/scripts/copytoImgPlus.m
3. Save mij.jar and imagej-matlab-0.7.2 in the ImageJ folder
4. Save OCTAVA.ijm into the folder ...ImageJ\Fiji\macros\OCTAVA Tools. You may need to create the folder named ¨OCTAVA Tools¨. 

*******************************************
Pre-compiled MATLAB app and standalone versions of OCTAVA are available at https://sourceforge.net/projects/octava/


********************************************
Rights and Permisisons

OCTAVA employs the following existing software packages

1. ImageJ (pulblicly available) 
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
