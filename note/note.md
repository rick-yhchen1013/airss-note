# Usage 
The airss calculation can be triggered by the following command:
airss.pl -keep -press 1 -max 60 -mpinp 16 -seed C2

The C2 is seed name of input file. the example can be found at the
"examples" directory.

See examples/2.1 for more detail


# airss call flow
airss.pl -> buildcell -> castep_relax -> cellsym -> run CASTEP -> castep2res

# The `ca -r` command

# How airss works with external libraries
