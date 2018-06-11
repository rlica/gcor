# gcor
Automatic Gain Correction relative to a reference spectra based on minimised chi sqared (for LaBr3(Ce) detectors)
Details about the code were published in [UPB Sci. Bull. A 78 (2016)](https://www.scientificbulletin.upb.ro/rev_docs_arhiva/rezce9_838206.pdf)

## Author
Razvan Lica (IFIN-HH), razvan.lica@nipne.ro

Contributors: C. Costache and N. Marginean (IFIN-HH)

====


## Requirements:
 1. The graphical interface (x11) for gnuplot   
 2. The Gnu Scientific Library (GSL)   
 
 * (Ubuntu) `sudo apt-get install gnuplot-x11 libgsl-dev`    
 * (osx)    `sudo port install gnuplot gsl`     
 
## Compiling: 
 Use the Makefile (thanks to N. Marginean)
```make gnu```
  
## Usage:
 1. Place the binary calibrated spectra (using GMATCH if needed) from all the runs and all the detectors (GASPware projections) in the same folder
 2. Run `gcor`.
 3. Modify the default values from the `gcor.settings` file to suit your needs (highly recommended).
 4. Press "Enter" to cycle each run, or type "a" and press "Enter" so that the code will cycle automatically all the files without using the graphical display.
 5. The correction coefficients are written in the `gcor.cal` file.

## Other tools:

### GMATCH (Gain Matching)
Gain Matching between a calibration run and run taken during the experiment because
the gain changes linearly due to different rates

This program will require the following arguments:
  - the calibration coefficients extracted from the calibration runs (.cal / .mcal)
  - a list ( [Det#] [Energy] [Channel] ) from the spectra during the experiment (only one peak for each detector)
It will output the calibration coefficients as a .cal file. 

*Usage:          ```gmatch LE-152Eu.cal list.dat```


### PINT (Polynomial Intersection)
Used for GASPWare .mcal files: estimates the coordinate of intersection between calibration polynomials
in order to establish the correct corresponding fit regions.

*Usage:          ```pint4 file.mcal```

