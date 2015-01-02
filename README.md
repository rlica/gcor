gcor
====

Automatic Gain Correction relative to a reference spectra based on minimised chi sqared

You must install:
 1. The graphical interface (x11) for gnuplot
 2. The Gnu Scientific Library (GSL)
 
 (Ubuntu) sudo apt-get install gnuplot-x11 gsl-dev
 (osx)    sudo port install gnuplot gsl
 
 
 
 Compile with: gcc -o gcor4 gcor4.c -lm -lgsl -lgslcblas (linux)
               gcc -o gcor4_osx gcor4.c -lm -L/opt/local/lib -lgsl -lgslcblas (osx)
