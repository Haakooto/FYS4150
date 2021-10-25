# Project3 - Numerical Study of the Penning Trap with N Particles
Sara Pernille Jensen, Alec Elías Sigurðarson, Håkon Olav Torvik

This repository contains our code for Project 3 of FYS4150 fall 2021.

## Running code
The c++ code can be compiled and run though the `makefile`, by using the command "`$ make all`" in the terminal. This uses c++ version 20, so might not work if this is not install. Changing this to version 14 or 17 should also correctly compile the code. By calling this, all the implemented experiments will be run. This can take a while. Specific experiments can be run by calling it in the main-function. The results are saved in the directory `outputs`, so make sure this is present.

The plots are made by calling "`$ python plot.py`". If all experiments have been run, all plots in the report should be recreated.

## Code structure
Two classes are in the code, PenningTrap and Particle. The codes for both classes are stored in the folder src/ in separate files. The header file for both classes is in include/PenningTrap.hpp.

The main code for running all the experiments is in main.cpp, which returns output in the folder outputs/. These results can then be plotted using the python script.