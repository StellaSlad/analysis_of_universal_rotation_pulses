# Analysis-of-broadband-Universal-Rotation-NMR-pulses

Scripts and functions for the analysis of broadband Universal Rotation pulses used in NMR spectroscopy. 

You can create graphs showing the Transfer Function of the pulses, or the total rotation axis for different offsets. 
The script will analyze all pulses in the current directory, and 

A sample input file showing our pulse format is available. 
It is important not to store any other file types in the same directory!

Planned improvements:

- add function that can read .bruker input files
- add a function in the beginning that automatically recognizes files with other endings and puts them into the subdirectory
"other_files" or creates output message "Detected a file format that cannot be processed, file name: ... . Please move this file to a different
directory".
