#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <sys/io.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sched.h>
#include <fcntl.h>
#include <complex.h>
#define PI 3.1415926536
#define VELC 2.99792458e08
#define NDAT 2048

// simplified for EDGES3 - just VNA calibration using SOL on 8-T switch
void cabsparms(complex double,complex double, complex double, complex double,complex double, complex double, complex double *, complex double *, complex double *);
complex double agilent(double,double,double);


int main(int argc, char *argv[])
{
  complex double s11[NDAT],s22[NDAT],s1221[NDAT];                     // s-parms for lab
  complex double Tfopen[NDAT],Tfshort[NDAT],Tfload[NDAT],Tfant[NDAT];  // data from field
  complex double ss11ant[NDAT];
  complex double s11open[NDAT],s11short[NDAT],s11load[NDAT];  // S11 data for Agilent SOL
  complex double T;
  double freqq[NDAT],freq,db,ph,res,loadps,openps,shortps;
  int i,j,k,m,n,csv;
  char name[255],buf[255];
  FILE *file3,*file;
  m = 0;  n = 1; T = 0; loadps = 38; openps = shortps = 33; res = 50; csv = 0;
  for(i=0;i<argc;i++){
  m = 0;
  sscanf(argv[i], "%79s", buf);
  if (strstr(buf, "-Tfopen")){ sscanf(argv[i+1],"%s",name); m=1;}
  if (strstr(buf, "-Tfshort")){ sscanf(argv[i+1],"%s",name);m=2;}
  if (strstr(buf, "-Tfload")){ sscanf(argv[i+1],"%s",name);m=3;}
  if (strstr(buf, "-Tfant")){ sscanf(argv[i+1],"%s",name);m=4;}
  if (strstr(buf, "-res")) sscanf(argv[i+1], "%lf",&res);
  if (strstr(buf, "-loadps")) sscanf(argv[i+1], "%lf",&loadps);
  if (strstr(buf, "-openps")) sscanf(argv[i+1], "%lf",&openps);
  if (strstr(buf, "-shortps")) sscanf(argv[i+1], "%lf",&shortps);
  if (strstr(buf, "-csv")) csv = 1;
  if(m){
  if ((file3 = fopen(name, "r")) == NULL) {
      printf("%s error here\n",name);
      return 0;
       }
   j = k = 0;
   while (fgets(buf, 255, file3) != 0) {
   if(strstr(buf,"DB")) k = 2;
   if(k || csv){
   freq = 0;
   if(csv) sscanf(buf, "%lf,%lf,%lf",&freq,&db,&ph);
   else sscanf(buf, "%lf %lf %lf",&freq,&db,&ph);
   if(freq > 0) {
   if(freq > 1e6) freqq[j] = freq/1e6;
   else freqq[j] = freq;
   T = pow(10,db*0.05)*cexp(I*ph*PI/180.0);
   if(csv) T = db + ph*I;
   if(m==1) Tfopen[j] = T;
   if(m==2) Tfshort[j] = T;
   if(m==3) Tfload[j] = T;
   if(m==4) Tfant[j] = T;
//  printf("m %d freq %f T %f %f %f\n",m,freq,creal(T),cimag(T),cabs(T));
   j++;
   n = j;
   }
   }
   }
   fclose(file3);
   }
   }
    sprintf(name,"s11.csv");
    if ((file = fopen(name, "w")) == NULL) {
        printf("cannot open %s:\n",name);
        return 0;
    }
   fprintf(file,"BEGIN\n");
   for(k=0;k<n;k++){

   if(freqq[k] > 1e6) {s11open[k] = agilent(freqq[k]*1e-6,1e9,openps); s11short[k] = agilent(freqq[k]*1e-6,0.0,shortps); s11load[k] = agilent(freqq[k]*1e-6,res,loadps);}
   else               {s11open[k] = agilent(freqq[k],1e9,openps); s11short[k] = agilent(freqq[k],0.0,shortps); s11load[k] = agilent(freqq[k],res,loadps);}

    cabsparms(s11open[k],s11short[k],s11load[k],Tfopen[k],Tfshort[k],Tfload[k],&s11[k],&s1221[k],&s22[k]); 
    ss11ant[k] = (Tfant[k] - s11[k])/(s1221[k]+s22[k]*(Tfant[k] - s11[k])); // get lna s1
    // SGM: changed write format from %f to %1.16e to increase precision.
    fprintf(file,"%1.16e,%1.16e,%1.16e\n",freqq[k],creal(ss11ant[k]),cimag(ss11ant[k]));
   }
   fprintf(file,"END\n");   
   fclose(file);
   return 0;
}


void cabsparms(complex double Topen, complex double Tshort, complex double Tload, 
    complex double Tobsopen,complex double Tobsshort, complex double Tobsload, complex double *s11, complex double *s1221, complex double *s22)
{
  complex double a00,a01,a02,a10,a11,a12,a20,a21,a22;
  complex double aa00,aa01,aa02,aa10,aa11,aa12,aa20,aa21,aa22;
  complex double b0,b1,b2,d;

   b0 = Tobsopen;
   b1 = Tobsshort;
   b2 = Tobsload;

   a00 = 1;  a01 = Topen;  a02 = Topen * Tobsopen;
   a10 = 1;  a11 = Tshort; a12 = Tshort * Tobsshort;
   a20 = 1;  a21 = Tload;  a22 = Tload * Tobsload;

    d = a00*a11*a22 + a10*a21*a02 + a20*a01*a12 - a20*a11*a02 - a10*a01*a22 - a21*a12*a00;

    aa00 =  (a11*a22-a21*a12)/d;
    aa01 = -(a01*a22-a21*a02)/d;
    aa02 =  (a01*a12-a11*a02)/d;
    aa10 = -(a10*a22-a20*a12)/d;
    aa11 =  (a00*a22-a20*a02)/d;
    aa12 = -(a00*a12-a10*a02)/d;
    aa20 =  (a10*a21-a20*a11)/d;
    aa21 = -(a00*a21-a20*a01)/d;
    aa22 =  (a00*a11-a10*a01)/d;

    *s11      = aa00*b0 + aa01*b1 + aa02*b2;
    *s1221    = aa10*b0 + aa11*b1 + aa12*b2;
    *s22      = aa20*b0 + aa21*b1 + aa22*b2;
    *s1221 += (*s11)*(*s22);
}

complex double agilent(double freq,double res,double delayps)
{ double delay,loss,rload;
  complex double Zcab,gl,Z,T,Tload;
// Agilent approx. follow
    loss = 2.30*1e9; 
//    printf("loss %f Gohm/s\n",loss/1e9);
    Zcab = 50 + (1 - I)*(loss/(2*2*PI*freq*1e6))*sqrt(freq*1e6/1e9);
//    delay = 33e-12;
    delay = delayps*1e-12;
    rload = res;
//    printf("Zcab %f %f\n",creal(Zcab),cimag(Zcab));
    gl = loss*(delay/(2*50))*sqrt(freq*1e6/1e9) + I*(2*PI*freq*1e6*delay + (loss*delay/(2*50))*sqrt(freq*1e6/1e9));
    T = (rload - Zcab)/(rload + Zcab);
    T = T*cexp(-2*gl);
    Z = Zcab*(1+T)/(1-T);
    Tload = (Z-50)/(Z+50);
//  printf("Tload %f %f\n",creal(Tload),cimag(Tload));
    return Tload;
}


