# NuclearPhysicsCodes
This README is seperated into subsections based on each of the files in this project. Each section describes the file and what changes you can make and what lines would need changed for that change to be processed.

AnimatedGraph.py
-------------------
Takes in five h5 files and plots their abundance versus mass number as an animation through time/Temp
Note: This plotting code uses the time values to align all the plots, so trajectories that start later will also start changing on the plot later.
Changing this to Temp in the future would probably be more helpful

To change the h5 files lines 22-26 set the paths to the h5 files 
Note: It is assumed all the h5 files are named WinNet_data.h5

Line 173 is used to set the xlim and ylim paramaters with a tuple
E.g. xlim=[0, 95] means x axis will go from 0 -> 95
     ylim=[1e-30, 1] means y axis will go from 1e-30 -> 1
Line 174 sets the scale of the plot options are 'log' and 'linear

After the animation has been shown on screen, you may close it and it will prompt you if you want to save the result.
This will save the result as a mp4 to your current working directory.


AnimatedNuclideTable.py
-------------------
Takes in a h5 file that contains tracked nuclei Winnet data and plots the isotopes and their abundances on the nuclide table in the form of an animation.

There are five spots to load h5 files but only one is used at a time. Line 60 can be modified to select a different one of the five written in h5 files.
Alternatively, you can also change the path or filename to the h5 data on lines 12-17

Lines 57 and 58 set the start and end points of the x and y axis respectively in the form of a tuple.
E.g. xr = (0,150) means x axis will go from 0 -> 150
     yr = (0,90) means y axis will go from 0 -> 9


After the animation has been shown on screen, you may close it and it will prompt you if you want to save the result.
This will save the result as a mp4 to your current working directory.


FlowGraph.py
-------------------


MainPlotter.py
-------------------


ManualGraphing.py
-------------------


ReverseEngineering,py
-------------------


StatFlowGraph.py
-------------------


TabulaterRateCalculator.py
-------------------


Track_Nuclei.py
-------------------


TrahUnitConv.py
-------------------


TrajectoryCalculator
-------------------


calcRates.oy


vo-proc.par
