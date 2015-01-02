#include<stdio.h>
#include <stdlib.h>
#include<math.h>
#include <gsl/gsl_multifit.h>

int chNum, detNum, runref, runstop;
int sweep, degree=1, maxChiSq;  
int left[100], right[100], chan[100], width[100];
char name[20];
int regions;
double sens;
int low, high;

typedef struct Data2Fit { 
  int ch;
  int chShift;
  double err;
} Data2Fit;


int    initialize() {
  FILE *settings;
  int i=0;
  
  if (fopen("gcor_settings.txt", "rt")) 
  {
    settings = fopen("gcor_settings.txt", "rt");
    printf("Settings taken from 'gcor_settings.txt'.\n");
    
  }
  else
  {
    printf("Printing default settings in 'gcor_settings.txt'. Modify them as you please :) \n");
    if(fopen("gcor_settings.txt", "wt"))
    {
      settings = fopen("gcor_settings.txt", "wt");
      
    }
    else { printf("ERROR: Cannot create 'gcor_settings.txt'\n"); exit(0); }
    fprintf(settings, \
"Channels=8192\n\
DetNum=11\n\
RefFile=L0.0001\n\
EndFile=L0.0200\n\
Sweep=150\n\
Sensitivity=4.0\n\
MaxChiSq=10\n\
LEFT | RIGHT\n\
  30     100\n\
 110     180\n\
 190     310\n\
 320     430\n\
 440     570\n\
 580     780\n\
 790     990\n");
    fclose(settings);
  }
  settings = fopen("gcor_settings.txt", "rt");
  fscanf(settings, \
"Channels=%d\n\
DetNum=%d\n\
RefFile=%5[^.].%04d\n\
EndFile=%5[^.].%04d\n\
Sweep=%d\n\
Sensitivity=%lf\n\
MaxChiSq=%d\n\
LEFT | RIGHT\n", &chNum, &detNum, name, &runref, name, &runstop, &sweep, &sens, &maxChiSq);
  
  // %5[^.] means read at most five characters or until a dot is encountered. google is the best :)
  
  
  printf("------------\nChan\t%d\n\
DetNum\t%d\n\
RefFile\t%2s.%04d\n\
EndFile\t%2s.%04d\n\
Sweep\t%d\n\
Sensitivity=%lf\n\
MaxChiSq\t%d\n\
LEFT | RIGHT\n\
------------\n", chNum, detNum, name, runref, name, runstop, sweep, sens, maxChiSq);


  
  while(fscanf(settings, "%d %d", &left[i], &right[i])!=EOF) 
  {
    if (left[i]>right[i] || left[i] < 0 || right[i] > chNum) 
    {
      printf("ERROR! gcor_settings.txt: Region limits %d\t%d are wrong\n", left[i], right[i]);
      exit(0);      
    }
    printf("%4d%8d\n", left[i], right[i]);
    width[i]=right[i]-left[i];
    chan[i]=(right[i]+left[i])/2;
    i++;
  }
  
  return i; //the number of regions
}

void   centerMax(double *refSpec) { //will center each region at the highest channel inside it
  
  int i, j;
  int max=0, jmax=0;
  int shift;
  for (i=0; i<regions; i++)
  {
    max=0, jmax=0;
    for (j=left[i]; j<right[i]; j++)
    {
      if (refSpec[j]>max) 
      {
	max=refSpec[j];
        jmax=j;
      }
    }
    chan[i]=jmax;
    width[i]=right[i]-left[i];
    shift = (sweep*chan[i])/chan[regions-1];
    if (chan[i]-width[i]/2<low+shift) chan[i]+=(low+shift)-(chan[i]-width[i]/2); //to avoid going beyond the low limit
    if (chan[i]+width[i]/2>high-shift) chan[i]-=(width[i]/2+chan[i])-(high-shift); //to avoid going above the high limit
  }
  for (i=1; i<regions; i++) 
  {
    if (chan[i]-chan[i-1]<width[i]/2) width[i]=0; //to remove overlapping regions
  }
  
}

int    readData(int **rawSpec, int *irun, int *irunref) {
  
  int i, runCount=0;
  FILE *fi;
  char infile[20];
   
  
  for(i=1; i<=runstop; i++)
  { 
  
  sprintf(infile, "%2s.%04d", name, i);
  if(!fopen(infile, "rb")) 
  {
    if (i==runref) {printf("ERROR: Cannot open reference %2s.%04d\n", name, runref); exit(0); }
    else continue;
  }
  else fi=fopen(infile, "rb");
  fread(rawSpec[runCount], sizeof(int),detNum*chNum,fi);
  irun[runCount]=i;
  if (i==runref) *irunref=runCount;
  runCount++;
  
  fclose(fi);
  }
  
  if (runCount<2) {printf("ERROR: Could not read more than the reference run\n"); exit(0);}
  else printf("----> %d runs were successfully read\n", runCount);
  
  return runCount;
        
       
}

void   polynomialfit(int n, int degree, double *xi, double *yi, double *ei, double *chisq, double *coeff) {

  //Modified by R. Lica. Taken from
  //        http://rosettacode.org/wiki/Polynomial_regression
  //        http://www.gnu.org/software/gsl/manual/html_node/Fitting-Examples.html

  int i, j;
  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;

  X = gsl_matrix_alloc (n, degree);
  y = gsl_vector_alloc (n);
  w = gsl_vector_alloc (n);

  c = gsl_vector_alloc (degree);
  cov = gsl_matrix_alloc (degree, degree);

  for (i = 0; i < n; i++)
    {
      
      
      gsl_matrix_set (X, i, 0, 1.0);
      for(j=0; j < degree; j++)
      {
	gsl_matrix_set(X, i, j, pow(xi[i], j));
	
      }
      
      
      gsl_vector_set (y, i, yi[i]);
      gsl_vector_set (w, i, 1.0/(ei[i]*ei[i]));
    }

  {
    gsl_multifit_linear_workspace * work 
      = gsl_multifit_linear_alloc (n, degree);
    gsl_multifit_wlinear (X, w, y, c, cov,
                          chisq, work);
    gsl_multifit_linear_free (work);
  }
/*
#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

  {
    printf ("# best fit: Y = %g + %g X + %g X^2\n", 
            C(0), C(1), C(2));

    printf ("# covariance matrix:\n");
    printf ("[ %+.5e, %+.5e, %+.5e  \n",
               COV(0,0), COV(0,1), COV(0,2));
    printf ("  %+.5e, %+.5e, %+.5e  \n", 
               COV(1,0), COV(1,1), COV(1,2));
    printf ("  %+.5e, %+.5e, %+.5e ]\n", 
               COV(2,0), COV(2,1), COV(2,2));
    printf ("# chisq = %g\n", chisq);
  }
*/

  for (i=0; i<degree; i++) coeff[i]=gsl_vector_get(c,i);
  
  gsl_matrix_free (X);
  gsl_vector_free (y);
  gsl_vector_free (w);
  gsl_vector_free (c);
  gsl_matrix_free (cov);

}

double normalize(double *spec1, double *spec2) {
  
  
  int i;
  double sum1=0, sum2=0;
  double norm=0;
  for (i=low; i<high; i++) { sum1+=spec1[i]; sum2+=spec2[i]; }
  if (sum1>sum2)
    for (i=low; i<high; i++)
    {
      norm = sum1/(sum2+1);
      spec2[i]=spec2[i]*norm;
    }
  else
    for (i=low; i<high; i++)
    {
      norm = sum2/(sum1+1);
      spec1[i]=spec1[i]*norm;
    }
  //printf("Norm. = %.3lf\t", norm);
  return norm;
}

void   smooth(double *spec1, double *spec2, int minW, int maxW) {
  
  int i, j;
  int ilow=low+minW/2, ihigh=high-maxW/2;
  double sum1=0, sum2=0, width=0;
  
  for (i=ilow; i<ihigh; i++)
  { 
    sum1=0, sum2=0, width=0;
    for (j=-(minW+(maxW-minW)*i/(ihigh-ilow))/2; j<(minW+(maxW-minW)*i/(ihigh-ilow))/2; j++)
    {
      sum1+=spec1[i+j];
      sum2+=spec2[i+j];
      width++;
    }
    spec1[i]=sum1/width;
    spec2[i]=sum2/width;
  
    //fprintf(fo, "%d\t%d\n", i, spec1[i]); 
  }
  
  
  
}

void   deriv(double *spec1, double *spec2, int der_width) {
  
  
  int i;
  double der1[high-der_width], der2[high-der_width];
  for (i=low; i<high-der_width; i++) {
    der1[i]=(spec1[i+der_width]-spec1[i])/der_width;
    der2[i]=(spec2[i+der_width]-spec2[i])/der_width;
  }
  
  for (i=low; i<high-der_width; i++) {
    //fprintf(out, "%d\t%lf\t%lf\t%lf\t%lf\n", i, spec1[i], der1[i], spec2[i], der2[i]);
    
    spec1[i]=der1[i];
    spec2[i]=der2[i];
     
  }
  
  
  
}

double chisq(double *spec1, double *spec2, int l1, int l2, int n) { //chisq of two regions
  
  int i;        
  double sum=0;
  double c2=0;
  
  for (i=0; i<n; i++)
  {
    c2 += pow(( spec1[l1+i]- spec2[l2+i]), 2);
    sum+= abs(spec1[l1+i])+abs(spec2[l2+i]);
  }
  
  return c2/sum;  
}

int    area(double *spec, int l, int n) {
  int i;
  double a=0;
  for (i=l; i<l+n; i++) a+=spec[i];
  return a;
}
    
double ratio(double *spec, int l, int n) { //checks if there is a significant change in the region; 
  int i;
  double min=spec[l], max=spec[l];
  double ratio=0;
  double limit;
  for (i=l; i<l+n; i++)
  {
    if (spec[i] < min ) min = spec[i];
    if (spec[i] > max ) max = spec[i];
  }
  
  
  ratio = max - min;
  limit = sqrt(max)+sens*sqrt(min);
  if (ratio > limit ) return ratio;
  else return 0;
    
}

void   shift(double *smoothSpec, double *spec1, double *spec2, struct Data2Fit *shData) { // finds shift and error for each region by sweeping and minimizing chiSq
 
  int i, j;
  double minchi;
  double signif;
  int max_shift;  
  for(i=0; i<regions; i++) if(width[i]>0)
  {
    shData[i].ch=0;
    shData[i].chShift=0;
    shData[i].err=0;
    
    minchi=chisq(spec1, spec2, chan[i]-width[i]/2, chan[i]-width[i]/2, width[i]);   
    signif=ratio(smoothSpec, chan[i]-width[i]/2, width[i]);
    
    if (signif)
    {
      max_shift = 5 + (sweep*chan[i])/chan[regions-1]; //version 1.4 
      
      for (j=-max_shift; j<max_shift; j++) 
	if(chan[i]-width[i]/2+j > low && chan[i]+width[i]/2+j < high)
      {
	if( minchi >= chisq(spec1, spec2, chan[i]-width[i]/2, chan[i]-width[i]/2+j, width[i]))
	{
	  
	  minchi=chisq(spec1, spec2, chan[i]-width[i]/2, chan[i]-width[i]/2+j, width[i]);
	  
	  shData[i].ch = chan[i]; 
	  shData[i].chShift = -j;
	  shData[i].err = (15+sqrt(abs(j)))/sqrt(signif);
	  
	  //shData[i].err = 2;
	  //shData[i].err = sqrt(minchi)*area(spec2, low, high)/(area(spec2, i*SHIFT, WIDTH[i]));
	  //printf("%d\t%d\t%f\t%d\n", shData[i].ch, shData[i].chShift, shData[i].err, WIDTH[i]); //testing
	  //if (i>0) if (abs(shData[i].chShift - shData[i-1].chShift) > 50) shData[i].ch = -1;
	
	}
	
      }
      
	
    }
    else shData[i].ch = -1;
    
  }
  
}

int    performFit(struct Data2Fit *shData, double *chisq, double *coeff) {
  
  int i, j=0, n=0;
  for (i=0; i<regions; i++)
    if(shData[i].ch > 0) n++;
  double xi[n], yi[n], ei[n];
  
  for (i=0; i<regions; i++) 
    if(shData[i].ch > 0)
    { 
      xi[j]=shData[i].ch;
      yi[j]=shData[i].chShift; 
      ei[j]=shData[i].err;      
      j++;
      
      //printf("%d\t%d\t%f\n", shData[i].ch, shData[i].chShift, shData[i].err); //testing
    
      
    }
  if (j<2) n=0; //not enough points to perform fit
  else polynomialfit(n, degree, xi, yi, ei, chisq, coeff);
  
  /*
  for(i=0; i < degree; i++) {
    printf("%lf\n", coeff[i]);
  }
  printf("%g\n", *chisq/n); // normalised chi sqared = chisq/degrees of freedom
  */
  return n;
    
}

void   gnuplot(FILE * gnuplotPipe, int irun, int idet, struct Data2Fit *shData, int n, double chisq, double *coeff,double norm) {
  
  
  
  char title[200], plot[200], labels[200], legend[200];
  sprintf(title, "set title \"DET#%02d RUN#%04d  Chisq=%.2lf  Norm=%.2lf\"", idet, irun, chisq/n, norm);
  if(degree == 2) sprintf(plot, "plot [0:1500][-%d:%d]'gcor_fit-data.temp' using 1:2:3 with yerrorbars, %lf + %lf*x", sweep, sweep, coeff[0], coeff[1]);
  else if(degree == 3) sprintf(plot, "plot 'gcor_fit-data.temp' using 1:2:3 with yerrorbars, %lf + %lf*x + %lf*x*x", coeff[0], coeff[1], coeff[2]);
  else if(degree == 4) sprintf(plot, "plot 'gcor_fit-data.temp' using 1:2:3 with yerrorbars, %lf + %lf*x + %lf*x*x %lf*x*x*x", coeff[0], coeff[1], coeff[2], coeff[3]);
  sprintf(labels, "set xlabel \"Channel\" \n set ylabel \"Shift\"");
  sprintf(legend, "set key left top box 3");
  
  FILE * temp = fopen("gcor_fit-data.temp", "wt");
  //FILE * gnuplotPipe = popen ("gnuplot", "w");
  int i;
  for (i=0; i < regions; i++)
    {
      if(shData[i].ch > 0)
	fprintf(temp, "%d %d %lf \n",  shData[i].ch, shData[i].chShift, shData[i].err); //Write the data to a temporary file
	fflush(temp);
    }
   

    fprintf(gnuplotPipe, "%s \n %s \n %s \n %s \n",labels, legend, title, plot); //Send commands to gnuplot one by one.
    fflush(gnuplotPipe);
  fclose(temp);
    
    
}

void   gnuplot_spec(FILE * gnuplotPipe, double *spec1, double *spec2, char *title, int num) {
  
  
  int i;
  FILE *out;
  char filename[200], labels[200], plot[200];
  sprintf(filename, "gcor_spec%d.temp", num);
  sprintf(labels, "set xlabel \"Channel\" \n set ylabel \"Counts\"");
  sprintf(plot, "plot '%s' u 1:2:3 with lines linecolor variable", filename);
  
  out=fopen(filename, "wt");
  
  for (i=low; i<high-30; i++) fprintf(out, "%d\t%lf\t1\n", i, spec1[i]);
  fprintf(out, "\n");
  for (i=low; i<high-30; i++) fprintf(out, "%d\t%lf\t2\n", i, spec2[i]);
  
  fflush(out);
  fclose(out);
  
   

  fprintf(gnuplotPipe, "%s \n unset key \n set title \" %s \" \n %s \n ", labels, title, plot); //Send commands to gnuplot one by one.
  fflush(gnuplotPipe);
     
}

void   writeCal(FILE *fo, double *coeff, int irun, int idet) {
 
  
  
  int deg;
  fprintf(fo, "%5d%5d%5d", irun, idet, degree);
  for (deg=0; deg<degree; deg++)
      {
	if (deg==0) fprintf(fo, "%9.3f", coeff[deg]);
	else if (deg==1) fprintf(fo, "%10.6f", coeff[deg]);
	else fprintf(fo, "%15.6E", coeff[deg]);
      }
  fprintf(fo, "\n");
  fflush(fo);
    
	
}

double testFit(double *refSpec, double *Spec) {
 
  int i;
  double *refTemp, *Temp, *smoothSpec, chi=0, signif;
  refTemp = calloc(chNum, sizeof(double));
  Temp = calloc(chNum, sizeof(double));
  smoothSpec = calloc(chNum, sizeof(double));
  for (i=low; i<high; i++) { refTemp[i]=refSpec[i]; Temp[i]=Spec[i]; }
  
  smooth(refTemp, Temp, 3, 15);
  for (i=low; i<high; i++) smoothSpec[i]=Temp[i];
  normalize(refTemp, Temp);
  //deriv(refTemp, Temp, 3);
  
  chi = chisq (refTemp, Temp, low+100, low+100, high-200);
    
  
  return chi; 
  
}

double calibrate(double *refSpec, double *Spec, double *coeff) {
 
  
  
   int i, j;
   double calj;
   double *calSpec;
   double chi;
   calSpec = calloc(chNum, sizeof(double));
   
   for (j=low; j<high; j++)
   {
     calj=0;
     for (i=0; i<degree; i++)
     {
      calj += coeff[i] * pow(j, i);
     }
     
     calSpec[ (int)calj ] += (1-abs((int)calj - calj))*Spec[j];
     calSpec[(int)calj+1] += abs((int)calj - calj)*Spec[j];
          
   }
   for (j=low+1; j<high; j++)           //to remove the gaps
   {
     if (calSpec[j]<0.5*calSpec[j-1]) calSpec[j]=(calSpec[j-1]+calSpec[j+1])/2;
     if (calSpec[j]>1.5*calSpec[j-1]) calSpec[j]=(calSpec[j-1]+calSpec[j+1])/2;
     Spec[j]=calSpec[j];
   }
   
   chi=testFit(refSpec, calSpec);
      
   return chi;  
}



