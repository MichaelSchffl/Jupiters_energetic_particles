# Jupiters energetic particles

 -- This is a short summary and explanation file for this project -- 

## Introduction
In this project we want to know more about energetic particles in Jupiter's vast magnetosphere. Specifically, we want to learn about their energy distribution with respect to the magnetic equator. We also want to learn about the angle they form with the local magnetic field, called the pitch angle (PA). In order to achieve this, we use data from the JEDI instrument of the Juno spacecraft (https://link.springer.com/article/10.1007/s11214-013-0025-3). This instrument detects high energetic ions in its orbit around Jupiter. Every 53 days it performs one full orbit. The data for this project were obtained from the JMIDL software of the JHUAPL and can be downloaded from this repository.
Ultimately, our goal is to speculate about a heating mechainsm that accelerates ions in Jupiter's middle magnetosphere by looking at the calculted distribution functions.


### About the data
The twelve orbits (Orbit01 to Orbit12) contained in the Orbits folder in this repository contain timespans where data have been measured and were preprocessed to .d2s files with another software. This code is specified to work with the .d2s files provided. The full calibrated and uncalibrated data from 2016 through 2019 is publicly available here: https://pds-ppi.igpp.ucla.edu/mission/JUNO/JNO/JEDI.


### How to run the code
For this project, only the Jupyter Notebook is needed and no additional software. The Jupyter Notebook can be installed here [https://jupyter.org/install] or with the Anaconda distribution (https://www.anaconda.com/products/individual).
The main script is commented so that it should be clear what is the purpose of each programming step. This is a quick walkthrough how to run the code and how to interprete the results it is providing:
- Download the repository and navigate to the folder you store it in the Jupyter environment. Run the main.ipynb script in Jupyter notebook.
- A cell can be run with ctrl-enter or shift-enter to jump to next cell
- When asked for the Orbit time and ion species, enter the orbit number as dd-dd_mm_, so e.g. for Orbit09 that is 17-20_10_ and the ion species as either LoTOFxE, HiTOFxE or HiResIon depending on which ion dataset is used in the orbit dataset you have chosen. 
- The functions.py file contains all functions written to support the code. These functions are imported at the beginning of the script. Additionally, numpy, scipy, matplotlib and datetime, ipywidgets and future are libraries that are needed to run the code 

### Tasks
 - The data sets provide particle intensities as a function of their energy (keV) and as a function of the pitch angle, along with the time the measurements were taken. As a first step, the orbit and ion data sets must be chosen and imported as well as the magnetic field data. The magnetic field data give information about the current spacecraft location w.r.t the magnetic equator.
 - These intensities need to be allocated to the times in the magnetic field data so that we know, where exactly the spacecraft is located. For that, the year, month, day, hour, minute and second columns have to be rounded to the nearest integer and compared to every value in the magnetic field data set. Since the magnetic field data are much higher resolved, meaning there is more rows than in the energy or pitch angle data, these data sets must be interpolated to the same length as the magnetic field data set.
 
 
 ### Results
 The results shown here are far from being exhaustive of what can be shown with data analysis done here, but give a good overview of what the data mean and how a data interpretation can look like.<br/>
 The following image shows the dipole latitude of Jupiter in blue and the magnetic field pattern in the same plot in green. Notice that the magnetic field values are uncalibrated.<br/>
<img src="https://github.com/MichaelSchffl/jupiters_energetic_particles/blob/master/images/Orbit09_magnetic_field.png"> <br/>
Following the code, the imported density values can be plotted for each species, along with the magnetic field measured by the spacecraft.<br/>
<img src="https://github.com/MichaelSchffl/jupiters_energetic_particles/blob/master/images/Orbit09_particle_densities.png"> <br/>
For the changepoint search to see where the spacecraft crosses the current sheet, the residual minimization method yields the following resul.t<br/>
<img src="https://github.com/MichaelSchffl/jupiters_energetic_particles/blob/master/images/Orbit09_residual_minimization.png"> <br/>
