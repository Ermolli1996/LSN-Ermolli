=================================================
| Numerical Simulation Laboratory - Exercise 04 |
=================================================
|                    Author:                    |
|                  Marco Ermolli                |
|        marco.ermolli@studenti.unimi.it       	|
|        Universita' degli Studi di Milano      |
=================================================

Legend of files:

Code is insterted the file named "MolDyn_NVE.cpp" .

Compile the code with terminal command (" g++ -o MolDyn_NVe.exe MolDyn_Nve.cpp ")

Executable is labelled with ".exe" extension.

Execute with "./< name_executable >" command.

Clean all "output*" files with "./clean.sh" command.

Input file are called "input.dat", and config.0. There are also other input files examples for solid, liquid and gas for input.dat and a config.fcc for config.0 .

Outputs are: "output*", "config.final" and "old.final".

Folder named "Equi" include output file of a solid equilibration (Ex 4.1) copied from the main folder.

Folders named "Solid", "Liquid" and "Gas" include all the results of the real simulation (Ex 4.3) copied from the main folder.

Folder named "Frames" is empty: it will become full of outputs file of molecular coordinates if you include in the main function che commented line ( "ConfXYZ(nconf)" )
