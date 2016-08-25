
This software implements the shallow-water equations with spherical harmonics.

On a long-term basis, it is intended to merge this development with SWEET.


Short guide for utilization:

First, install SHTNS library for Spherical Harmonics transformations:
	$ cd ./local_software
	$ ./install_shtns.sh
	$ cd ..


Setup environment (Don't forget the '.' at the very beginning):
	$ . ./local_software/env_vars.sh
	OR
	$ source ./local_software/env_vars.sh


Compile the software:
	$ make release
