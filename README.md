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
This is the main plotting code for plotting Abudances and Flows as a function of temperature. It plots four datasets at each timestep
It works by taking in a list of temperature values and finds the abundances and reaction flows at that temperature. 
And instead of plotting as an animation as in the previous codes, it outputs each temp step as its own png file.

To change the h5 files lines 28-33 set the paths to the h5 files 
Note: It is assumed all the h5 files are named WinNet_data.h5

Line 784 is used to set the temp values to evaluate at in units of Gk
E.g. list(np.around(np.linspace(7, .1, 100),3)) will start at 7GK and go to .1 Gk in 100 even steps, the around function is used to round the temp values to 3 decimal places

Lines 788 and 789 set the start and end points for the x and y axis respectively. The first entry being the starting point and the second the ending.

Lines 906 through 909 set the Titles for the four plots. 

Line 939 is where you set the path for all of the images to be placed in. 

If you prefer a video format for the plots, the naming system has been setup so you can use ffmpeg to combine all of the images as frames into a video.

At the very end of Line 939 the dpi= paramter sets the clarity of the image. The higher the value the more detailed the output image (enabling you to zoom in and such) 
I found 300 to be a good value where it doesn't take too long to save but it still fairly clear even zoomed in. 


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


calcRates.py
-------------------


vo-proc.par
-------------------
