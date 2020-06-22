=================================================
| Numerical Simulation Laboratory - Exercise 07 |
=================================================
|                    Author:                    |
|                  Marco Ermolli                |
|        marco.ermolli@studenti.unimi.it       	|
|        Universita' degli Studi di Milano      |
=================================================

Legend of files:

Code is insterted the file named "Monte_CarloNVT.cpp" .

Compile codes with "make" command.

Executable is labelled with ".exe" extension.

Execute with "./< name_executable >" command.

Clean all outputs files with "./clean.sh" command.

Input file are called "input.dat", and "config.0". There are also other input files examples for solid, liquid and gas for "input.dat" and a config.fcc for "config.0" .

Outputs are: "output*", "config.final" and "seed.out".

Folder named "Solid", "Liquid" and "Gas" include output files of all the simulations in the three phases copied from the main folder.

Folder named "MolDyn_NVE" includes codes and results of Molecular Dynamic simulations. 

Folder named "Frames" is empty: it will become full of outputs file of molecular coordinates if you include in the main function che commented line ( "ConfXYZ(nconf)" )

All other files are used for generating ramdon numbers.
