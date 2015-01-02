/***************************************************************************************************************
*                         GainCor for LaBr3:Ce spectra - R. Lica - Dec. 2013                                   *
*  Automatic Gain Correction relative to a reference spectra based on minimised chi sqared (not a peak finder) *
****************************************************************************************************************
 
 You must install:
 1. The graphical interface (x11) for gnuplot
 2. The Gnu Scientific Library (GSL)
 
 (Ubuntu) sudo apt-get install gnuplot-x11 gsl-dev
 (osx)    sudo port install gnuplot gsl
 
 
 
 Compile with: gcc -o gcor4 gcor4.c -lm -lgsl -lgslcblas (linux)
               gcc -o gcor4_osx gcor4.c -lm -L/opt/local/lib -lgsl -lgslcblas (osx)
 
 -----------------------------
Version 1.4 (6 Feb 2014)
-----------------------------
-the maximum shift (sweep) changes linearly with the channel (not a fixed value any more)
-fixed some mistakes related to the low and high limits 
-----------------------------
Version 1.3 (20 Dec 2013)
-----------------------------
-can perform the algorithm twice in order to correct strongly shifted spectra.
-removed higher degree polynomial calibration.
-removed the derivative spectra display. To enable, uncomment lines 180-181.
-added the display of recalibrated spectra
-every region will be auto-centered on the highest channel inside it
-if the fit is not good (normalised chisq greater than a default value), the program will 
  warn the user and skip the output to 'gcor.cal'.

-----------------------------
Version 1.2 (14 Dec 2013): 
-----------------------------
-the user will set which regions to consider -> autosearch is not reliable yet
   Disadvantage: the same regions will be considered for all the detectors 
                 -> it is recommended that the reference run to have all 
                    the detectors aligned
-added plotting of raw and derivative spectra
-any run can be set as reference
-----------------------------
Version 1.1 ( 5 Dec 2013): 
-----------------------------
-variable width, smoothing of spectra, derivative.
-----------------------------
Version 1.0 (29 Nov 2013): 
-----------------------------
-can read one detector from several runs, normalization of spectra.

*****************************************************************************************************************/
  


#include "gcor4_definitions.h"

int main(int argc, char **argv) {

  
  int i, j, k;
  printf("\n======= GCOR v1.4 - Automatic Gain Correction =======\n\tR. Lica, IFIN-HH, February 2014 \n\n");
  
  
  //Initializing 
    
  regions=initialize();
  degree++; // in the program, 'degree' represents the number of coefficients of the polynomial, NOT the degree of the polynomial
  //if (degree<2 || degree>5) {printf("ERROR: Degree must be between 1 and 4"); exit(0);}
  low=chan[0]-width[0]/2;
  high=chan[regions-1]+width[regions-1]/2+sweep;
  
  
  struct Data2Fit shData[regions];
  int onetime=0, cal_again=0;
  double refSpec[chNum], Spec[chNum], smoothSpec[high-low];
  int nData=0;
  double coeff[degree], coeff_old[degree], chisq, norm;
  char answer;
  char title[200];
  FILE * gnuplotPipe  = popen ("gnuplot", "w");
  FILE * gnuplotPipe2 = popen ("gnuplot", "w");
  FILE * gnuplotPipe3 = popen ("gnuplot", "w");
  FILE * gnuplotPipe4 = popen ("gnuplot", "w");
  FILE *fo;
  char outfile[20];
  sprintf(outfile, "gcor.cal");
  fo=fopen(outfile, "wt");
  
  //Reading the data
  int irun[runstop], irunref, idet, runCount;
  int **rawSpec; // same as  'int rawSpec[runstop][detNum*chNum];' but can have any size
  rawSpec = calloc(runstop, sizeof(int *));
  for(i = 0; i < runstop; i++) rawSpec[i] = calloc(detNum*chNum, sizeof(int));
  runCount = readData(rawSpec, irun, &irunref); 
  
  double chiTest, chiTest_cal;
  
  
  //GCOR main code
  
  for(idet=0; idet<detNum; idet++)
  {
    onetime=0;
    printf("\n---------------------\nDetector #%02d\n---------------------\n", idet);
      
  
  for(i=0; i<runCount; i++) {
    if(irun[i]!=runref)
  {
    cal_again=0; 
  
    for (j=0; j<chNum; j++)
    {
      refSpec[j]=rawSpec[irunref][idet*chNum+j];
      Spec[j]=rawSpec[i][idet*chNum+j];
    }
    
    chiTest=testFit(refSpec, Spec);
  
    if (answer!='a') 
    {
    sprintf(title, "Raw Spectra. ChiSq = %.2lf", chiTest);
    gnuplot_spec(gnuplotPipe3, refSpec, Spec, title, 1); 
    }
  
    if (onetime==0) { centerMax(refSpec); onetime++;}
  
  again:
  for (j=0; j<regions; j++) 
    {
      shData[j].ch=0;
      shData[j].chShift=0;
      shData[j].err=0;
    }
  smooth(refSpec, Spec, 10, 30);
  for (j=low; j<high; j++) smoothSpec[j]=Spec[j];
   
  norm = normalize(refSpec, Spec);
  deriv(refSpec, Spec, 3);
  //sprintf(title, "Derivative");                               // Uncomment this lines in order to display 
  //gnuplot_spec(gnuplotPipe2, refSpec, Spec, title, 3);           //  the derivative spectra
  
  shift(smoothSpec, refSpec, Spec, shData); 
  nData = performFit(shData, &chisq, coeff);
  if (nData == 0) 
  {
    printf("\n%2s#%02d.%04d: No output in gcor.cal: Could not extract data suitable for fit.\t", name, idet, irun[i]);
    //fprintf(fo, "%5d%5d%5d%9.3f%10.6f\n", irun[i], idet, 2, 0.0, 1.0);
    if (answer!='a')
    {
      printf("Going to next? [y]/n/a\t");
      answer = getchar();
    }
    goto skip;
    
  }
  
  if (answer!='a' && cal_again==0) gnuplot(gnuplotPipe, irun[i], idet, shData, nData, chisq, coeff, norm);
  coeff[1]+=1;        
  
  if (cal_again==0)   //performing again the algorithm
  {
    coeff_old[0]=coeff[0];
    coeff_old[1]=coeff[1];
    for (j=0; j<chNum; j++)
    {
      refSpec[j]=rawSpec[irunref][idet*chNum+j];
      Spec[j]=rawSpec[i][idet*chNum+j];
    }
    calibrate(refSpec, Spec, coeff_old);
    cal_again++;
    goto again;
  }
  else                //updating the coefficients
  {
    coeff[0]=coeff[0]+coeff[1]*coeff_old[0];
    coeff[1]=coeff_old[1]*coeff[1];
  }
  
  for (j=0; j<chNum; j++)
    {
      refSpec[j]=rawSpec[irunref][idet*chNum+j];
      Spec[j]=rawSpec[i][idet*chNum+j];
    }
  chiTest_cal=calibrate(refSpec, Spec, coeff);
  
  if (answer!='a') //disable the display in automatic mode
  {
    
    sprintf(title, "Recalibrated Spectra. ChiSq = %.2lf", chiTest_cal);
    gnuplot_spec(gnuplotPipe4, refSpec, Spec, title, 2);
    }
  
  
  
  if (chiTest_cal > chiTest && chiTest_cal<maxChiSq) 
  {
    fprintf(fo, "%5d%5d%5d%9.3f%10.6f\n", irun[i], idet, 2, 0.0, 1.0);
    if (answer!='a')
    {
      printf("\n%2s#%02d.%04d: Recalibration is not required (cal/raw = %.2lf/%.2lf) \t", name, idet, irun[i], chiTest_cal, chiTest);
      printf("Going to next? [y]/n/a\t");
      answer = getchar();
    }
  }
  
  else if (chiTest_cal>maxChiSq && chiTest>maxChiSq && chiTest/chiTest_cal<2)
  {
    printf("\n%2s#%02d.%04d: No output in gcor.cal: Could not perform a good fit (cal/raw = %.2lf/%.2lf) \t", name, idet, irun[i], chiTest_cal, chiTest);
    if (answer!='a')
    {
      printf("Going to next? [y]/n/a\t");
      answer = getchar();
    }
    
  }
  
  else 
  { 
    writeCal(fo, coeff, irun[i], idet);
    if (chiTest_cal>maxChiSq && chiTest>maxChiSq)
    {
      printf("\n%2s#%02d.%04d: Needs improvement (cal/raw = %.2lf/%.2lf) \t", name, idet, irun[i], chiTest_cal, chiTest);
      if (answer!='a')
      {
	printf("Going to next? [y]/n/a\t");
        answer = getchar();
      }
    }
    else if (answer!='a') 
      {
	printf("\n%2s#%02d.%04d\t OK (cal/raw = %.2lf/%.2lf)\t", name, idet, irun[i], chiTest_cal, chiTest);
        printf("Going to next? [y]/n/a\t");
        answer = getchar();
      }    
  }
  
  skip:
  if (answer=='n') exit(0);
  
    
  }
  
  else fprintf(fo, "%5d%5d%5d%9.3f%10.6f\n", runref, idet, 2, 0.0, 1.0); //each detector from the reference run is set as reference
  }
  
    
  }
   
  exit(0); 
  
  
}