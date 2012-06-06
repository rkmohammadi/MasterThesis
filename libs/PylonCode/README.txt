This archive contains the MATLAB code for the inference in the pylon segmentation model described in:

V. Lempitsky, A. Vedaldi, and A. Zisserman
A PYLON MODEL FOR SEMANTIC SEGMENTATION.
Advances in Neural Information Processing Systems (NIPS), Granada, 2011.

//////////////////////////

QUICK START:

Download the QPBO code (see "Dependencies" below). Then type "pylonExample" in the MATLAB prompt.

//////////////////////////

CONTENTS:

"pylonInference*Class.m" - implements inference for the 1 class, 2 class and multiclass case, as described in the paper. See also the comments in the beginning of each function for the usage suggestion.

"pylonSetup.m" - compiles the QPBO code into a mex file (required to run the inference).

"pylonExample.m" - demonstrates the usage of the pylon model. It takes the example image (from the Berkeley Software Dataset) along with the segmentation tree (computed with the Arbelaez et al. algorithm) and the brushes for four image parts and segments the image into the four parts using a) a multiclass pylon model, b) a multiclass pylon model without pairwise terms, c) a flat CRF code. Note that the majority of the time is spent on visualization. The time spent on the single call to the multiclass pylon inference will be displayed.

"pylonVisualize.m" - a rather inefficient code for the visualization of the pylon segmentation result.

//////////////////////////

DEPENDENCIES:

1) The pylon inference calls the QPBO algorithm inside it. Therefore, you need to download the QPBO code (version 1.3) from Vladimir Kolmogorov's website (http://pub.ist.ac.at/~vnk/software.html). Unpack the code files (*.cpp,*.h,*.inc) into the main folder.
Note, that although this version of the code calls the QPBO code, it is possible to use the submodular graph cuts after a simple variable substitution (described in the paper). Since this substitution exists, QPBO is guaranteed to label all nodes (as long as the pairwise terms are submodular).

2) To run the example provided with the code, you need to install the VLFeat toolbox (www.vlfeat.org). This code was tested with the VLFeat version 0.9.14.

3) To run the example provided with the code, you need to have the MATLAB Statistics and Image Processing toolboxes installed.
