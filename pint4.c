//  P(olynomial)INT(ersect) v1.3 - R. Lica, IFIN-HH, Feb 2013
//  Used for .mcal files -> estimates the coordinate of intersection between calibration polynomials
// in order to establish the correct corresponding fit regions

//Updates:
//24 Jan 2014 (v1.4) - Warning if no intercept is found (Limited to 10.000.000 iterations)
//12 Dec 2013 (v1.3) - improved output of .mcal file
//11 Dec 2013 (v1.2) - if intercept is < 0 the program will return 0 as the 'newlimit' in order to discard the first polynomial
// 5 Feb 2013 (v1.1) - correct approximation of intersecton (int = double + 0.5)


#include<stdio.h>
#include <stdlib.h>
#include<math.h>
#define EPS 0.00001
#define RANGE 10
#define MAXDEG 5



double F(double x, int deg, double *pol1, double *pol2)
  {
    double ans=0; 
    int i;
    for (i=0; i<deg; i++) ans+=(pol1[i]-pol2[i])*pow(x,i);
    return ans;
    
  }

void readMaxDeg(FILE *fi, int *maxdeg)   //v1.3
{
  int line=0, col=0, regions=0, i, j, k;
  int deg[30][10]; //max 30 lines(detectors), 10 columns(regions)
    for (i=0; i<30; i++) for (j=0; j<10; j++) deg[i][j]=0;
  
  while (fscanf (fi, "%*d %*d %d %d", &regions, &deg[line][col])!=EOF) 
  {
    if( maxdeg[col] < deg[line][col] ) maxdeg[col] = deg[line][col];
    
    for (k=0; k<deg[line][col]; k++) fscanf (fi, "%*f");	  
    fscanf (fi, "%*d"); // limit1);
	
	
    for (i=0; i<regions-1; i++)
    {
      col++;
      fscanf (fi, "%d", &deg[line][col]);
      for (k=0; k<deg[line][col]; k++) fscanf (fi, "%*f");
      fscanf (fi, "%*d"); //limit2);
      if( maxdeg[col] < deg[line][col] ) maxdeg[col] = deg[line][col];
    }
    line++;
    col=0;    
  }
  
  
}
  
int main(int argc, char **argv)
{

  int i, j, k, step=0, newlim;
  int runnum, detnum, regions, deg1, deg2, deg, limit1, limit2;
  double pol1[MAXDEG], pol2[MAXDEG], x1, x2, x3, f1, f2, t;
  char outfile[120];
  FILE *fi;
  FILE *fo;
  
  
  if (argc < 2) {
    printf("\nPINT v1.4 (Jan 2014) - ERROR:  .mcal file required as argument:\n\n");
    exit(0);
  }
  

      
  fi=fopen(argv[1], "r");
  sprintf(outfile, "pint_%s", argv[1]);
  fo=fopen(outfile, "wt");
  
  /// version 1.3
  int reg=0;
  int maxdeg[10]; //max 10 columns(regions)
  for (j=0; j<10; j++) maxdeg[j]=0;
  readMaxDeg(fi, maxdeg);
  fclose(fi);
  /// 
    
  for (j=0; j<MAXDEG; j++) { pol1[j]=0; pol2[j]=0; }
  
  fi=fopen(argv[1], "r");
  while (fscanf (fi, "%d %d %d %d", &runnum, &detnum, &regions, &deg1)==4) {
    
	
	reg=0;
	fprintf(fo, "%5d%5d%3d%5d", runnum, detnum, regions, deg1);
	for (k=0; k<deg1; k++) {
          fscanf (fi, "%lf", &pol1[k]);
	  if (k==0) fprintf(fo, "%9.3f", pol1[k]);
	  else if (k==1) fprintf(fo, "%10.6f", pol1[k]);
	  else fprintf(fo, "%15.6E", pol1[k]);
	  if (k==deg1-1) 
	  {
	    while(k+1<maxdeg[reg]) { fprintf(fo, "               "); k++;}      //v1.3
	    break;
	  }
	  
	}
	  
	fscanf (fi, "%d", &limit1);
	reg++;
	
	
    for (i=0; i<regions-1; i++) {
	  fscanf (fi, "%d", &deg2);
	  for (k=0; k<deg2; k++) fscanf (fi, "%lf", &pol2[k]);
	  fscanf (fi, "%d", &limit2);
	  
	  x1=limit1-RANGE;
	  x2=limit1+RANGE;
	  
	  if (deg1>deg2) deg=deg1;
	  else deg=deg2;
	  
	  int control=0;
      do {       // Secant method
        f1=F(x1, deg, pol1, pol2);
        f2=F(x2, deg, pol1, pol2);
        x3=x2-((f2*(x2-x1))/(f2-f1));
        x1=x2;
        x2=x3;
        if(f2<0)    t=fabs(f2);
        else    t=f2;
	control++;   //version 1.4
	if (control>10000000) {printf("Warning! Detector #%02d Regions %d-%d: Functions do not intersect.\n", detnum, i, i+1); break;}
        } while(t>EPS);
	  
	  
	  newlim=(x3+0.5);
	  
	  ///////////// v1.2 update
	  if (newlim<0) 
	  { 
	    newlim=0;
	    printf("\nWarning! Detector #%02d: First limit is %0.2lf. Returned 0 in the .mcal file", detnum, x3);	    
	  }
	  /////////////
	  
	  fprintf(fo, "%6d%5d", newlim, deg2); //x3 is the new limit
	  
	  limit1=limit2;
	  deg1=deg2;
	  for (j=0; j<MAXDEG; j++) pol1[j]=0;
	  for (k=0; k<deg2; k++) {
	    pol1[k]=pol2[k];
	    if (k==0) fprintf(fo, "%9.3f", pol2[k]);
	    else if (k==1) fprintf(fo, "%10.6f", pol2[k]);
	    else fprintf(fo, "%15.6E", pol2[k]);
	    if (k==deg2-1)
	    {
	      while(k+1<maxdeg[reg]) { fprintf(fo, "               "); k++; }  //v1.3
	      break;
	    }
	    
	  }
	  
	  for (j=0; j<MAXDEG; j++) pol2[j]=0;
	  step++;
	  reg++;
	  
	  }
	fprintf(fo, "%6d\n", limit2);
	}  
	  
	printf("\nPINT v1.4 (Jan 2014) - %d updated limits written in: '%s'\n\n", step, outfile);
    exit(0);
}
