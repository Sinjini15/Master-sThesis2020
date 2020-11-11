# Master-sThesis2020


This file explains in detail how to run the code by combining files from both the MATLAB code folder and the Main Code Folder.

Please save both folders in the same location

Main Code:

- Contains files scripted in R and contains the main file "bradycardia_detection.R"
- Contains supporting functions that are used in the main file.

Instructions to run bradycarida_detection.R:

The code asks for user input initially. Make sure to input the correct parameters.
The parameters asked from the user are:
- sampling frequency of the signals (can be obtained from the dataset)
- split of the training/validation/testing set
- value of probability of false alarm
- the type of kernel to be used

Follow the comments in the .R and run the relevant MATLAB functions where directed. Subsequent sections may depend on the MATLAB file being run.
Save the paths to the outputs of the .mat files as they are called later

All graphs can be plotted by the plot commands contained at the bottom of the main code.

MATLAB Code:

- Contains 4 files of the type cross_val_<insert kernel name> that computes the leave one out cross validation for each type of kernel
- Depending on what kernel is to be used, run the corresponding .m file. The file saves a .mat output which is called to R
- Also contains the pan-tompkin function that will be used to detect the peaks. This function is called when event_gen.m is run and is an existing library. 



