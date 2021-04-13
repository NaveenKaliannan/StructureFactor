
#!/bin/bash

make
./exe traj.xyz 100 1 15.6404 15.6404 15.6404 384 rdffile


#1) executable
#2) trajectory file (only xyz format accepted)
#3) trajectory length
#4) time step in fms
#5) box length in x direction
#6) box length in y direction
#7) box length in z direction
#8) number of atoms in the whole box
#9) output file (prints rdf)







