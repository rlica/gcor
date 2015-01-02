# gcor
Automatic Gain Correction relative to a reference spectra based on minimised chi sqared (for LaBr3(Ce) detectors)

## Author
Razvan Lica, IFIN-HH, razvan.lica@nipne.ro

====


## Requirements:
 1. The graphical interface (x11) for gnuplot   
 2. The Gnu Scientific Library (GSL)   
 
 * (Ubuntu) `sudo apt-get install gnuplot-x11 gsl-dev`    
 * (osx)    `sudo port install gnuplot gsl`     
 
## Compiling: 
 * (linux) `gcc -o gcor4 gcor4.c -lm -lgsl -lgslcblas`
 * (osx)   `gcc -o gcor4_osx gcor4.c -lm -L/opt/local/lib -lgsl -lgslcblas`

## Usage:
 1. Place the binary spectra from all the runs and all the detectors (GASPware projections) in the same folder
 2. Run `gcor`.
 3. Modify the default values from the `gcor.settings` file to suit your needs.
 4. Press "Enter" to cycle each run, or type "a" and press "Enter" so that the code will cycle automatically all the files without using the graphical display.
 5. The correction coefficients are written in the `gcor.cal` file.
