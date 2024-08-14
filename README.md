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
This is the main plotter for graphs showing trajectory charachtersitics or other large scale calculation values. It can be used to plot Abudance versus Mass number entropy vs time and more.
It grabs the data from 5 h5 Files the paths for which are set on lines 10-14

Line 9 also sets the number of output graphs in rows and columns so plots = [3, 3] would mean 3 rows 3 columns so 9 plots.
When using plots like this the axis object that is created is actually a numpy array of axis representing each of the plots so to axis any individual plot you have to do axis[y][x] 
where y is row number starting from the top at zero and x is the column number starting from the left at zero

The majority of the data is read straight out of the h5 files. To find the different children and subgroups of any particular group in a h5 file you can simply call .keys() which will return a key object that contains all the groups under the one you've opened.
Also when gathering the data from the file like in lines 551 through 555 for each of the five trajectories we plot the entropy verses Yn. f1out['entr'][:] takes the first h5 file looks under the 'entr' tag that is short for entropy and the [:] means to collect all values. In some cases you will also see me use [()]. They are equivalent. 

ManualGraphing.py
-------------------
The manual graphing file is set up at the moment to manually plot the flows of particular important reactions as a function of temperature.

You can change the paths for these plots on lines 12-16. Change the isotopes on lines 287-293. Extras could be added or removed.
And you can also change the indiviudal reactions watched on lines 313 through 319.
Plots can be edited on lines 358 and 364.

The code also repeats itself once to plot the second trajectory but that isn't strictly necesseary. 

ReverseEngineering,py
-------------------
This was the code I was trying to use to reverse engineer reaction rates. It builds a matrix of reactions rates and attempts to find an RREF to solve the system. After solving the system it will plot its calculated rates and flows versus those that we know to be true. It currently does not function very well

The input files h5 and Summed Flows File can be set on lines 32-46. The plot row and column numbers are set on line 696 and the xrange and yrange vales for plotting are set on 697 and 698.
SummedFlows input filename is set on line 704

This file is very experimental. At the moment it is doing a lot of unnecessary stuff for the purpose of debugging

StatFlowGraph.py
-------------------
This is very similar to the other flow graph file. Except it is meant for static flows it won't plot multiple plots or animations just a single output image that contains the releavant flows. It is currently setup to plot the total summed flows of each reaction in the calculation.

The xrange and yrange values for plotting are set on lines 754 and 755. Plot size set on line 764. SummedFlowsFilenames are set on line 838, 841, 845, and 848 for each of the four plotted trajectories. and the path to these files are set on line 32
The titles for each of the four plots are set on lines 860 through 863.

The output is set to both save to a file and show the user the output file path and name is set on line 892

TabulaterRateCalculator.py
-------------------
This file is used to take in important reactions find their tabulated rates from the reaclib website and automatically download format and save them to a file for use in Winnet. 

Line 9 is where you set the output path for the tabulated rates file. 
Line 182 is where each of the temperature values for the tabulated rates is set. These values are in GK, this line should not need to be edited it is set as the default for Winnet, but can be changed if necessary. 

The workflow for this file is meant to pretty simple for reaction rate you want to tweak call findRateIndex with a string of the reaction you want to find a tabulated rate for "cu59(p,a)ni56" for example is a valid string. Weak rates can be found by leaving the values inside the parenthesis empty "ni56(,)co56" for example would represent a Ni56 beta decay. After you have called this and grabbed all of its output values you can call findTabulatedRates which will actually return the rates. The rates are stored in a numpy array and you can use array or matrix operations to change the rates. As an example I have the rates being doubled on lines 188, 193, and 198 for each of the three reactions in the example. The line outString = appendRatestoOutfile(detName, vers, res, dir, Q, rates, outString) just finally adds all of the data to a string for future saving. Once you have done this for all the different reactions you want to change you can call cleanOutString and writeOutput on lines 202 and 203 which will format the string for printing and remove unnecessary data and sections. 


Track_Nuclei.py
-------------------
This file creates a list of isotopes for winnet to track in the track_nuclei paramater. 

The output path and filename are set on lines 3 and 4. Mass ranges are set ib kubes 6 and 7. There are three options for adding isotopes to be tracked. The first one is to grab all the elements in a certain mass range. The second is to grab all isotopes within the same mass range which works better and provides more isotopes, but at the extreme ends of the nuclide table either towards extreme neutron richness or extreme proton richness it will not count and add these isotopes to be tracked. Not normally an issue just something to note. The last option is to manually add the isotopes you want to track by their element number. These three methods are shown in lines 65, 66, 67 respectively. They can be commented in and out based on you're prefered method. 

TrahUnitConv.py
-------------------
This file was created because the Trajectories sent to me by Peli had the wrong units for the Temperature. It was in Kelvin while Winnet required them to be in GK. But since it has acquired more features including filtering out ttrajectories based on their charachterstics and plotting different trajectories together at the end so you can visualze their diffrenreces. It will also tell you how most of the trajectories failed incase you want to add or remove trajectories from your final collection.

The lines that set the path to the folder with all the trajectories and the output path are on lines 10 and 11 respectively. Lines 17 to 45 set the minimum and maximum values for all of the trajectory charachterstics. Note: The code only considers the initial value of these paramaters. So long term the values will diverge from the ones you used to reduce the sample. 

TrajectoryCalculator
-------------------
Follows the calculations in Otsuki et al to try and develop trajectory files. Currently experimental.

Input paramaters are lines 57 to 60.

Does not currently contain a feature for outputting traj files not consistent enough results yet.

calcRates.py
-------------------
This file calculates the total summed reaction rates for each trajectory after nucleosynthesis has been performed. Was orginally developed for use with the Reverse Engineering file. Warning: This file takes several minutes to run but does have sliders and features to help show you a time estimate.

Input h5 files are set on libes 14 through 24. Each trajecory is processed independently. The code first works through the entire dataset obtaining a complete list of reactions that take place and then goes through and sums them through time to try and determine the total flows generated. 

The output file name is ste on line 626 for trajecory 1 and similar for the rest. 

vp-proc.par
-------------------
The Winnet paramater file used for all of the nucleosyntheis calculations. 
