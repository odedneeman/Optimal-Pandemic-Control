# Code for "Integrating Macroeconomic and Public Health Impacts in Designing Social Planning Policies for Pandemic Response"

**Authors:** Ofer Cornfeld, Kaicheng Niu, Oded Neeman, Michael Roswell, Gabi Steinbach, Stephen J. Beckett, Yorai Wardi, Joshua S. Weitz, and Eran Yashiv
This repository contains the MATLAB code used to generate the results and plots for the epi-econ model described in our paper. The model evaluates the effects of various intervention policies during an infectious disease pandemic.

**Overview**

*	All results are saved as .mat files but can be fully reproduced using the provided MATLAB scripts.
*	The code for generating the figures is located in the main folder and can be run without any prerequisite computations.

**Folder Structure**

*	**Common**: Contains auxiliary functions for use in optimizations and simulations.
*	**ContinuousTime**: Contains code for continuous-time policy optimizations and general simulations.
*	**DiscreteTime**: Contains code for discrete-time policy optimizations and simulations.
*	**FeedbackControl**: Contains code for policy optimizations and simulations using a feedback control algorithm.
*	**accumHarm**: Conatins some results of welfare loss across diseases and policies.

