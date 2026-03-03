# Dynamics-Project---Analysis-of-a-6-rod-mechanism

Kinematic and dynamic analysis of a 6-rod mechanism and animate the mechanism using Matlab and Adams.

**Author:** Sajad Sarafzadeh
**Affiliation:** University of Tehran
**Contact:** sajad.sarafzadeh@gmail.com

---

## Project Overview
This repository contains the comprehensive kinematic and dynamic simulation of a 6-bar mechanism, developed as a Machine Dynamics course project. The mechanism is driven by a constant input angular velocity of 1 rpm on the driving link (rod 2). 

## Objectives & Tasks
* Calculate the Degree of Freedom (DOF) of the system, verifying that the system has 1 DOF.
* Formulate appropriate vector loops to mathematically define and extract the position, velocity, and acceleration of each linkage.
* Plot specific kinematic graphs, including the absolute velocity and acceleration of slider D, and the angular velocity and acceleration of rod CD.
* Perform a complete force analysis to determine the reaction forces acting on joints B, C, and D, as well as the required driving torque on rod 2.
* Develop an animation of the mechanism's continuous motion across all time steps.

## Methodology
* **Kinematic Analysis:** The geometric constraints were solved by establishing two closed vector loops. These position equations were then differentiated with respect to time to extract the velocity and acceleration profiles.
* **Dynamic Analysis:** The Newton-Euler equations of motion were derived for all individual rods and sliders. The mathematical model was constrained with a linear rod density of 1 kg/m, a 300 kg mass for the slider, and a kinematic friction coefficient of 0.25.

## Software & Tools Used
* **MATLAB:** Used as the primary environment for coding the mathematical models, numerically calculating the kinematics and dynamics at each time step, and outputting the mechanism's animation.
* **SAM:** Utilized to simulate the mechanism and provide a baseline for verifying the kinematic charts.
* **MSC ADAMS:** Employed to build a 3D model of the system to simulate and validate both the kinematic and dynamic (force/torque) calculations.

## Results & Conclusion
* The kinematic tracking (position, velocity, and acceleration curves) obtained from the custom MATLAB code demonstrated excellent alignment with the output from commercial software (SAM and ADAMS).
* The dynamic force and torque results showed slight discrepancies when compared to ADAMS. This is primarily because ADAMS analyzes the mechanism within a 3D environment, whereas perfectly modeling continuous friction and gravitational effects throughout the cycle in the custom MATLAB script introduces minor calculation variations.
