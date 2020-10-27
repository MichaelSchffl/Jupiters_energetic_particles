# Jupiters_energetic_particles

 -- This is a short summary and explanation file for this project -- 

## Introduction
In this project we want to know more about energetic particles in Jupiter's vast magnetosphere. Specifically, we want to learn about their energy distribution with respect to the magnetic equator. We also want to learn about the angle they form with the local magnetic field, called the pitch angle (PA). In order to achieve this, we use data from the JEDI instrument of the Juno spacecraft (https://link.springer.com/article/10.1007/s11214-013-0025-3). This instrument detects high energetic ions in its orbit around Jupiter. Every 53 days it performs one full orbit. The data for this project were obtained from the JMIDL software of the JHUAPL and can be downloaded from this repository.
Ultimately, our goal is to speculate about a heating mechainsm that accelerates ions in Jupiter's middle magnetosphere by looking at the calculted distribution functions.

### Task

 - The datasets provide particle intensities as a function of their energy (keV) and as a function of the pitch angle, along with the time the measurements were taken. As a first step, these intensities will be allocated to the times we know the spacecraft is ... 
 
 
 
 ### How to run the code
The main script is commented so that it should be clear what is the purpose of each programming step. This is a quick walkthrough how to run the code and how to interprete the results it is providing:
- download the repository and run the script main.ipynb on Jupyter notebook. I recommend installing Jupyter notebook with the Anaconda distribution.
- A cell can be run with ctrl-enter or shift-enter to jump to next cell
- When asked for the Orbit and ion species, enter the orbit number as dd-dd_mm_, so e.g. for Orbit09 that is 17-20_10_ and the ion species as either LoTOFxE, HiTOFxE or HiResIon depending on which ion dataset is used in the orbit dataset you have chosen. The twelve orbits (Orbit01 to Orbit12) contained in the Orbits folder in this repository contain timespans where data have been measured and were preprocessed to .d2s files with another software. This code is specified to work with the .d2s files provided. The full calibrated and uncalibrated data from 2016 through 2019 is publicly available here: https://pds-ppi.igpp.ucla.edu/mission/JUNO/JNO/JEDI [https://pds-ppi.igpp.ucla.edu/mission/JUNO/JNO/JEDI].
