
This software implements the shallow-water equations with spherical harmonics.

On a long-term basis, it is intended to merge this development with SWEET.


Short guide for utilization:

First, install SHTNS library for Spherical Harmonics transformations:
	$ cd ./local_software
	$ ./install_shtns.sh
	$ cd ..

Second, install lapack library for REXI:
	$ cd ./local_software
	$ ./install_lapack.sh
	$ cd ..

If you don't have FFTW3 installed, install it:
	$ cd ./local_software
	$ ./install_fftw3.sh
	$ cd ..


Setup environment (Don't forget the '.' at the very beginning):
	$ . ./local_software/env_vars.sh
	OR
	$ source ./local_software/env_vars.sh


Compile the software:
	$ make release



Run software e.g. with
	./build/sh_example T64 P0 N1 T0

Parameters:
	T64: spectral resolution 64
	P0:  with program P0: shallow-water equations
	N1:  non-linearities
	T0:  timestepping method 0 (RK4)


Plot results, e.g.:
	./plot.py prog_eta_*.csv


Make video:
	ffmpeg -framerate 25 -pattern_type glob -i 'prog_eta*.png' -c:v libx264  -pix_fmt yuv420p awesome_video_swe_eta.mp4
