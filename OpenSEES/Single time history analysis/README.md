This folder contains some of the files from the flexible rocking model developed in OpenSEES, for running a single time-history analysis.


The current .tcl file for running a single analysis is run_thisIsIt.tcl which contains some inputs to the overall system and runs the model file thisIsIt_cleanHouse.tcl.
The MATLAB file thisIsIt.m calls OpenSeesGNG.exe, which is the custom build of OpenSEES containing the GNG material model, and runs the tcl files mentioned above.
This MATLAB file also does post processing of the data. Some settings for the system may need to be carefully curated in both the MATLAB file and the tcl files.
