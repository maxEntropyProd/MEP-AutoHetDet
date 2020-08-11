This directory contains and example input file, example.inp.  Since model output can be many GB, the main output files, of which there are several, are not included here. However, the log file (which is the output to the consol during a run) as well as the summmary file are included.  This was run under v4.5, but the results should be the same under v4.7. The example was run on a cluster using 90 cores.

SLURM was use to handle batch runs, and the file example.slurm included here could be used to submit a job to a SLURM queue; however, the directories in the file likely need to be modified.

All output files are in ASCII, but with headers designed to be read in by Tecplot. A Tecplot file, example.lay, has also been included that would read all the outputs for plotting if you manage to get the program to run :)
