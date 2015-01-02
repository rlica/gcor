/***************************************************************************************************************
*                         GainMatch for LaBr3:Ce spectra - R. Lica - Feb. 2014                                   *
*                Gain Matching between a calibration spectra and spectra taken during experiment                 *
*                          (The gain changes linearly due to different rates)
****************************************************************************************************************
 
 *Compile with:   gcc -o gMatch gmatch.c -lm
 *This program will require the following arguments:
  - the calibration coefficients for the callibration spectra (.cal / .mcal)
  - a list ( [Det#] [Energy] [Channel] ) from the spectra during the experiment (only one peak for each detector)
 *It will output the calibration coefficients as a .cal file. 
 *Usage:          gmatch LE-152Eu.cal list.dat
            


*****************************************************************************************************************/

#include<stdio.h>
#include <stdlib.h>
#include<math.h>
#include <string.h> /* strrchr, strcmp*/

#define EPS 0.00001





double F(double x, int degree, double *coeff) {
    double ans=0; 
    int i;
    for (i=0; i<degree; i++) ans+=(coeff[i])*pow(x,i);
    return ans;
  }

double ChanFromEnergy (double *coeff, int degree, double energy) {
  
  
  int control=0;
  double f1, f2, x1, x2, x3, t;
  
  x1=(energy-coeff[0])/coeff[1]-100; //approx limits based on a linear dependence
  x2=(energy-coeff[0])/coeff[1]+100;
  
  do {       // Secant method
        f1=F(x1, degree, coeff)-energy;
        f2=F(x2, degree, coeff)-energy;
        x3=x2-((f2*(x2-x1))/(f2-f1));
        x1=x2;
        x2=x3;
        if(f2<0)    t=fabs(f2);
        else    t=f2;
	control++;   //version 1.4
	if (control>10000000) {printf("ERROR! Cannot find polynomial roots.\n"); break;}
        } while(t>EPS);
  
  return x3;
    
  
}

int readMCal(FILE *fi, double *coeff, int idet, double energy) {
  
  int degree, limit1, limit2;
  int regions=0, i, j, k, index;
  double chan;
  char string[1000];
  
  while (fscanf (fi, "%*d %d %d %d", &index, &regions, &degree)!=EOF){
    
    
    if (index == idet) {
      for (k=0; k<degree; k++) fscanf (fi, "%lf", &coeff[k]); 	  
      fscanf (fi, "%d", &limit1); 
      chan = ChanFromEnergy(coeff, degree, energy);
      if ( chan < limit1 && chan > 0 ) return degree;
      else for (i=0; i<regions-1; i++) {
        fscanf (fi, "%d", &degree);
        for (k=0; k<degree; k++) fscanf (fi, "%lf", &coeff[k]);
        fscanf (fi, "%d", &limit2);
        chan = ChanFromEnergy(coeff, degree, energy);
	if ( chan < limit2 && chan > limit1 ) return degree;
	else limit1=limit2;
      }
      
    }
    
    else fgets(string,sizeof string,fi); //reads a line from the file
  }
    
	
  
  
  printf("ERROR: Cannot find Detector %d\n", idet);
  return 0;
}

int readCal(FILE *fi, double *coeff, int idet, double energy) {
  
  int degree, index;
  int i, j, k;
  double chan;
  char string[1000];
  
  while (fscanf (fi, "%*d %d %d", &index, &degree)!=EOF){
    
    if (index == idet) {
      
      for (k=0; k<degree; k++) fscanf (fi, "%lf", &coeff[k]);	  
      return degree;
    }
    
    else fgets(string,sizeof string,fi); //reads a line from the file
  }
  
  printf("ERROR: Cannot find Detector %d\n", idet);
  return 0;
}

void writeGainCal(FILE *fo, double *coeff, int degree, int idet) {
 
  int deg;
  fprintf(fo, "%5d%5d%5d", 1, idet, degree);
  for (deg=0; deg<degree; deg++)
      {
	if (deg==0) fprintf(fo, "%9.3f", coeff[deg]);
	else if (deg==1) fprintf(fo, "%10.6f", coeff[deg]);
	else fprintf(fo, "%15.6E", coeff[deg]);
      }
  fprintf(fo, "\n");
  fflush(fo);
    
	
}


int main(int argc, char **argv) {
  
  
  int i, j, k, degree, index;
  double coeff[5], channel, list_channel, energy;
  char *ext, ext_mcal[]=".mcal", ext_cal[]=".cal";
  
  printf("\n---------\nGainMatch v1.0\n---------\n");
  if (argc < 3) {
    printf("Usage:  gmatch (.cal/.mcal file containing calibration coefficients) (file containing a list [Det#] [Energy] [Channel]) \n ");
    exit(0);
  }
  
  FILE *calFile, *listFile;
  
  if(fopen(argv[1], "rt")) calFile = fopen(argv[1], "rt");
  else { printf("ERROR: Cannot open %s\n", argv[1]); exit(0);}
  
  if(fopen(argv[2], "rt")) listFile = fopen(argv[2], "rt");
  else { printf("ERROR: Cannot open %s\n", argv[2]); exit(0);}
  
  ext = strrchr(argv[1], '.');    
  printf("Extension is %s\n", ext); // .cal or .mcal
  
  
  FILE *fo;
  char filename[100], calname[100];
  sscanf(argv[1], "%[^.]s", calname);
  sprintf(filename, "gMatch_%s.cal", calname);
  fo = fopen (filename, "wt");    
  printf ("Printing in %s\n\n", filename);
 
  
  while ( fscanf(listFile, "%d %lf %lf", &index, &energy, &list_channel)!=EOF ) { 
    
    
    printf("Det #%d\t", index);
    
    
    if (strcmp(ext,ext_mcal)==0) {  degree = readMCal(calFile, coeff, index, energy); }
    if (strcmp(ext,ext_cal)==0)  {  degree = readCal (calFile, coeff, index, energy);}
    
    channel = ChanFromEnergy(coeff, degree, energy);
    
    degree = 2;               // finding the gain change
    coeff[0]=0;
    coeff[1]=channel/list_channel;   
    printf("Gain = %lf %lf %lf\n", coeff[1], channel, list_channel);
    writeGainCal(fo, coeff, degree, index);
    
    rewind(calFile);
    
    
  }
  
  
  fclose(fo);
  fclose(calFile);
  fclose(listFile);
  
  
  exit(0);
 
  
  
}