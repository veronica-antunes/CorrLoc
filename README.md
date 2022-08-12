# CorrLoc

Antunes et al., (2022); Insights into the dynamics of the Nirano Mud Volcano through seismic characterization of drumbeat signals and V/H analysis,
Journal of Volcanology and Geothermal Research, 107619, ISSN 0377-0273, DOI: https://doi.org/10.1016/j.jvolgeores.2022.107619.

**Requirements:**
* Python3
* numpy version 1.19.5
* Obspy version 1.2.2
* scipy version 1.6.1
* matplotlib version 3.2.2

CorrLoc code was developed to locate weal emergent signals by using the waveform amplitudes instead of picks. The code backprojects the envelope of the cross-correlation function between station pairs to give a likelihood location solution. CorrLoc runs automatically through the continuous data and provides a location probability map for each window segment.

**Installattion:**

Install the package requirements. The code was tested with the versions provided above. If you face some package imcompatibilities, create a python environment and install the specific versions: https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/


**Usage:**

Import the different functions present in the scr/ folder and use them as described below. The code runs in 6 steps. Each step can be run independently for more flexibility. Each step is detailed explained in the publication ([link](https://www.sciencedirect.com/science/article/pii/S0377027322001500), see below).

```
import src.CorrLoc as CorrLoc
```


**Example:**

run_CorrLoc.py file exemplifies how to use the different functions and steps to locate a signal. It also includes an example on how to prepare the signal for location. We provide waveform examples (in wavs/) and the necessary input data (in loc_example/) to use the codes. The results are saved in loc_example/results. To run the example just open a terminal window and type:

```
python3 run_CorrLoc.py
```

To apply to real data just open the file and change the input parameters accordingly. The time window, filter and velocities might need to be tested and adapted for better solution result.


**Functions:**

```
grid(sta_file, x_start, x_end, y_start, y_end, step, outfile=None, fontsize=12, 
     verbose=True, path='')
```
First step: Creates the GRID and station coordinates
* sta_file - station file, with coordinates of the stations and the reference point for the (0,0) coordinate. Needs to be in the same format as loc_example/stations_file.txt
* x_start, x_end - start and end of the grid in the x direction, in meters
* y_start, y_end - start and end of the grid in the y direction, in meters
* step - step for grid calculations, in meters
* outfile - if a name is provided it will generate a map with the grid and the stations 
* fontsize - fontsize for plot lables, default is 12
* verbose - if True it will print messages to indicate the stage of the calculations. Default is True.
* path - if a path is inserted it will generate the files in the path. default is None (it will generate the files in the current directory)


```
crosscorr(wavs_file, delta_time_w, verbose=True, path=None)
```
Second step: Computes the cross-correlation envelope between station pairs for every time window defined (sliding window effect).
* wavs_file - waveform file (or path)
* delta_time_w - time segment window, in seconds
* verbose - if True it will print messages to indicate the stage of the calculations. Default is True.
* path - if a path is inserted it will generate the files in the path. default is None (it will generate the files in the current directory)


```
difftime(vel1, vel2, step, verbose=True, path=None)
```
Third step: Computes the diferential time for a list of velocities. It will calculate for all velocity values between vel1 and vel2 whithin the step defined. Allows testing different velocities in case the velocity is not known.
* vel1, vel2 - first (inclusive) and last (exclusive) velocity values for differential time calculations.
* step - step in velocity for differential time calculations.
* verbose - if True it will print messages to indicate the stage of the calculations. Default is True.
* path - if a path is inserted it will generate the files in the path. default is None (it will generate the files in the current directory)


```
envvalue(wavs_file, outfile=None, point=None, indv_norm=False, verbose=True, path=None)
```
Fourth step: Calculates the envelope value for every point of the grid, according to the respective diferential time.
* wavs_file - waveform file (or path)
* outfile - name of the outputfile in case we want to plot a specific point
* point - [int(i),int(j),int(p),int(w)], the given point to plot the CC envelope
			i,j is the grid point
			p is the station pair 
			w is the number of the time window
* indv_norm - if True, it will normalize each individual backprojections. Use in the case a station pair is dominating the location solution. Default is False.
* verbose - if True it will print messages to indicate the stage of the calculations. Default is True.
* path - if a path is inserted it will generate the files in the path. default is None (it will generate the files in the current directory)


```
stackall(final_norm=True, verbose=True, path=None)
```
Fifth step: Stacks all the individual back-projections into a final location solution.
* final_norm - if True, normalizes the final solution, for a normalized resulting values (from 0-1). Default is True.
* verbose - if True it will print messages to indicate the stage of the calculations. Default is True.
* path - if a path is inserted it will generate the files in the path. default is None (it will generate the files in the current directory)


```
plotall(name_folder, wavs_file, ntotal_stat=None, tr_norm=False, source_loc_file = None, 
        source_loc = None, out_format='.png', map_units = None, norm_tw=False, threshold=0.75, 
        show_loc_pos = True, final_loc_file='final_locations.list', verbose=True, path=None)
  ```
Sixth and last step: 	Plots the final likely location solution in a map and the waveforms.

* name folder - saves the figures in a certain location with a specific name
* wavs_file - waveform file (or path)
* ntotal_stat -to plot stations that are offline as an horizontal line with y=0, otherwise will just use the existing traces in the waveforms. Default is None.
* tr_norm - if True, normalizes the amplitudes of the waveform traces in the plot. Default is False.
* source_loc_file - file with different location positions (in case of a known source loc). The file must be in the same format as the station file, with the event name instead of the station name. Used to plot a specific position. Default is None.
* source_loc - Name of the location provided in the source_loc_file, to plot its position. This position will be plotted as a star and it can be used to test the code and the parameters in case of a known source.
* out_format - format of the output figures. 'png, 'pdf', etc. Default is 'png'.
* map_units - define the units of the map in the figure. 'm' for meters and 'km' for kilometers. If None, calculates it automatically considering 'm' for values below 1000 and 'km' for values above 1000 m. Default is None.
* norm_tw - normalizes the image function for each time window, for vizualization purposes. If set True it will norlmalize the amplitudes for each time window. USefull in presence of signals with different amplitudes. It will also normalize noise, so this parameter should be used with care. Default is False.
* threshold - plots the location solution only if the image function for the window is above the defined threshold. Default is 0.75 (assuming normalization in the previous step).
* show_loc_pos - plots the point of the grid with the maximum value (final location). Default True. If False it will not plot the point of the grid with the max amplitude value.
* loc_file - saves the final grid position and time ranges in a file. The final location is the the point of the grid with the maximum amplitude value. If more than one solution is acheaved, then it will plot all possible solutions.
* verbose - if True it will print messages to indicate the stage of the calculations. Default is True.
* path - if a path is inserted it will generate the files in the path. default is None (it will generate the files in the current directory)

Together with the functions to compute the locatation, there is also a set of functions to test each step. Each function is explained in detail in the  src/test_steps.py file. It is also provided an example of these functions in the example exercise.


**If you use this code please cite as:**

Verónica Antunes, Thomas Planès, Anne Obermann, Francesco Panzera, Sebastiano D'Amico, Adriano Mazzini, Alessandra Sciarra, Tullio Ricci, Matteo Lupi (2022); Insights into the dynamics of the Nirano Mud Volcano through seismic characterization of drumbeat signals and V/H analysis,
Journal of Volcanology and Geothermal Research, 107619, ISSN 0377-0273, DOI: https://doi.org/10.1016/j.jvolgeores.2022.107619. (https://www.sciencedirect.com/science/article/pii/S0377027322001500)

**Or Bibtex citation:**
```
@article{Antunes2022,
title = {Insights into the dynamics of the Nirano Mud Volcano through seismic characterization of drumbeat signals and V/H analysis},
journal = {Journal of Volcanology and Geothermal Research},
pages = {107619},
year = {2022},
issn = {0377-0273},
doi = {https://doi.org/10.1016/j.jvolgeores.2022.107619},
url = {https://www.sciencedirect.com/science/article/pii/S0377027322001500},
author = {Verónica Antunes and Thomas Planès and Anne Obermann and Francesco Panzera and Sebastiano D'Amico and Adriano Mazzini and Alessandra Sciarra and Tullio Ricci and Matteo Lupi}
}
```

