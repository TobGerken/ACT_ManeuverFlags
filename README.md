# ACT_ManeuverFlags

Author: Tobias Gerken, tug15@psu.edu

Version: 1.0 (2019_04_01)

This code is used to generate ManeuverFlag files for the [ACT-America](https://www-air.larc.nasa.gov/missions/ACT-America/) airborne research campaig. The data used corresponds to:

    ACT-America: L3 Merged In Situ Atmospheric Trace Gases and Flask
    Data, Eastern USA", https://doi.org/10.3334/ORNLDAAC/1593
    
Code is provided for documentation purposes only.

## Explanation of Maneuver Flags 

Maneuver Flags (MF) are assigned to give an indication of aircraft flight maneuvers during research flights. MFs are matched to the timestamp of the L3 dataset and take integer values between 0 and 6. 

| MF definition | MF Value |
|-|-|
|take off| 1| 
|in-line ascent| 2| 
|in-line descent| 3| 
|spiral up| 4| 
|spiral down| 5| 
|constant altitude legs| 6| 
|landing| 7| 

**Note: Maneuver Flags were assigned after the fact. Care should be taken when using MFs for data analysis to ensure that indentified maneuvers are sensible and that extracted data based on MFs is suitable for the analysis in question**

## Brief summary of Maneuver Flag generation procedure

In brief the algorithm assigns MFs for each flight-leg in the following sequence:

1. **Maneuver indentification not possible** if missing altitude information
2. **Airplane on ground** if altitude above ground level is below 25 meters
3. **Level leg identification**
     1. Level legs outside of typical boundary layer heights are identified based on dz/dt gradient threshold using altitude above sea level
     2. Near ground level legs are identified using altitude above ground level information
     3. Level legs of less than 5 minute duration are discarded
4. **Profile indentification** is accomplished by identfying local minima and maxima
     1. The aircraft heading is used to determine whether the plane undergoes at least 300 degrees of rotation during **spiral ascents and descents**
     2. Profiles are quality controlled to ensure that total altitude difference is more than 500m and that the mean gradient exceeds 1m/s 
5. **Quality control procedures** are applied given that lower-level legs can be noisy due to terrain changes. 
6. **Manual changes to MF** are made for some flight-legs and documented in the code for transparency/ reproducability 
7. **Quality flags** are assigned to indicate a first measure of confidence in the MF. QC = 0 indicates likely low quality of flag assignment. Note that QC=1 does not indicate good quality, but rather that initial non-exhaustive quality checks were passed. 

##  References/ Data Access:

Davis, K.J., M.D. Obland, B. Lin, T. Lauvaux, C. O'dell, B. Meadows, E.V. Browell, J.H. Crawford, J.P. DiGangi, C. Sweeney, .J. McGill, J. Dobler, J.D. Barrick, and A.R. Nehrir. 2018. ACT-America: L3 Merged In Situ Atmospheric Trace Gases and Flask Data, 
Eastern USA. ORNL DAAC, Oak Ridge, Tennessee, USA. [https://doi.org/10.3334/ORNLDAAC/1593](https://doi.org/10.3334/ORNLDAAC/1593)
 
