# CODEV

##General start-up steps:

1. After downloading the repository, open Matlab (R2011a or later; tested on R2014b).

2. Navigate to the 'codev' directory (this will ensure it is in your workspace).


##One-time installation steps:

1. In the command window, type 'init_mpt3' and follow the prompts to install MPT3. 
If there are any issues, consult: http://control.ee.ethz.ch/~mpt/3/Main/Installation


##To run the CODEV tool:

1. Make sure the desired input file is in the 'codev'-directory. For this example, we use the input file named 'inCruise.m'.

2. In the command window, type: '[reachtube,safeFlag,T] = verifyMPC(@inCruise);'

3. The results will be saved to the variables 'reachtube', 'safeFlag', and 'T'. 
If the system is robustly safe, 'safeFlag=1' and the 'reachtube' will have length(T) elements corresponding with the timestamps in 'T'. Otherwise, 'safeFlag=0' and 'reachtube' will be an empty object.

4. Reachtube results can be plotted as follows:
    * A single dimension (i) vs time: 'plotReachVsTime(reachtube,i,T)'
    * Two dimensions (i,j): 'plotReachVsTime(reachtube,[i,j],T)'
    * Two dimensions (i,j): 'plot(reachtube.projection([i,j]))'