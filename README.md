# Monitoring-Time-Dependent-Image-Processes-for-Detecting-Shifts-in-Pixel-Intensities

This repository contains the code for the simulation study conducted in the paper "Monitoring Time Dependent Image Processes for Detecting Shifts in Pixel Intensities" by Yarema Okhrin, Viktoriia Petruk, and Wolfgang Schmid.

## Abstract

The study addresses the challenge of detecting shifts in image processes, where pixel intensities are modeled using a spatial autoregressive process. Shifts manifest in changes in average intensities, and the primary goal is to detect these shifts as swiftly as possible post-occurrence. To effectively handle high-resolution images, the paper proposes a scalable technique focusing on the monitoring of regions of interest (ROIs). For shift detection, the study employs multivariate exponentially weighted moving average (EWMA) control schemes along with various control statistics. The effectiveness of this unique approach is showcased through an extensive simulation study. Recommendations are also provided for practitioners on selecting appropriate charts, setup, and calibration.

## Code Overview

This repository contains a photo used for the simulations and an R file that encompasses:

1. **Pre-Simulation Calculations**: Computational steps required prior to commencing the simulation study.
2. **Simulation Examples**: Demonstrations of simulations utilizing each type of control statistic described in the paper.

