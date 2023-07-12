MAKEFILE:

	make galaxysim: makes galaxysim.o to be used by the user


galaxysim.c: Simulates a galaxy collision of two galaxies of N particles where N is divisible by 10 for a user decided timestep, and simulation time.

By default it is set to utilize 8 threads, a user may edit the file header "int threads = 8" to change that, they may also edit the header to change the softener from 100pc or the accelcap from 4*10^-6.

desired input | double: timestep, int multiple of 10:number of particles, double: time to simulate |
output | nbodyleap-timestepnumberofparticlestimetosimulate.dat


The program utilizes a point-to-point method and is parallelized using openmp, the user must use a compiler like gcc or clang which has OpenMP support if they decide to compile it in another way, the MAKEFILE is recommended. 


plotit: generates an animation using the dat file in question with ffmpeg within Python
plottwo: saves each individual frame to /partials/ from a dat file for use with an external program like ffmpeg. It saves in the form "*.png"

plottwo is better for larger N or larger amounts of timesteps as you can stop anyywhere and it pick it up alter

plotit is better for shorter usecases as it is more flexible and allows a live rotatable view.