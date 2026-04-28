#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <sys/io.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <sched.h>
#include <fcntl.h>
#define  PI 3.1415926536  
#define TWOPI 6.28318530717958
#define NDATA 100000

complex double  cabl(double, double, complex double,complex double *, complex double *,double,double);

    int main(int argc, char *argv[])
{
    char fname[256];
    int j,k,nant,m;
    complex double Ta,s11,s12,s22,ss11[NDATA];
    double freq,fre,cablen,cabloss,cabdiel,re,im,freqq[NDATA];
    char buf[256],name[256];
    FILE *file2,*file3;
  nant = 1500;  cablen = cabloss = cabdiel = 0;
   if (argc > 1) sscanf(argv[1], "%s", name);
  for(m=0;m<argc;m++){
     sscanf(argv[m], "%79s", buf);
    if (strstr(buf, "-cablen")) {sscanf(argv[m+1], "%lf",&cablen);}  // length in inches
    if (strstr(buf, "-cabloss")) {sscanf(argv[m+1], "%lf",&cabloss);}  // loss correction in percent
    if (strstr(buf, "-cabdiel")) {sscanf(argv[m+1], "%lf",&cabdiel);}  // dielectric correction in percent
   } 
  if ((file3 = fopen(name, "r")) == NULL) {
      printf("%s error here\n",name);
      return 0;
       }
   j = k = 0;
   while (fgets(buf, 255, file3) != 0) {
   if(strstr(buf,"BEGIN")) k = 1;
   if(strstr(buf,"END")) k = 0;
   if(k>=2){      freq = 0;
   sscanf(buf, "%lf,%lf,%lf",&freq,&re,&im);
   if(freq){
   freqq[j] = freq;
   ss11[j] = re + im*I;
   j++;
   nant = j;
   }
   }
    k++;
   }
   fclose(file3);
    snprintf(fname, sizeof(fname), "c_%.253s", name);
    if ((file2 = fopen(fname, "w")) == NULL) {
        printf("cannot open file:%s\n", fname);
        return 0;
    }
    fprintf(file2,"BEGIN\n");
    for(j=0;j<nant;j++){
       freq = freqq[j];
       Ta=s11=s12=0;
       if(freq > 1e6) fre = freq/1e6; else fre = freq;
       cabl(fre,fabs(cablen)*2.54e-2/3e08,0,&s11,&s12,1+cabloss*0.01,1+cabdiel*0.01); // 0.15 ns cable in field cablen distance from ant input 
//       cabl(fre,cablen*2.54e-2/3e08,0,&s11,&s12,1+cabloss*0.01,1+cabdiel*0.01); // 0.15 ns cable in field cablen distance from ant input 
       s22 = s11;
       Ta = ss11[j];
       if(cablen > 0.0) Ta= s11+(s12*s12*Ta)/(1-s22*Ta);
       if(cablen < 0.0) Ta =(Ta-s11)/(s12*s12-s11*s22+s22*Ta);
//   Ta= s11+(s12*s12*Ta)/(1-s22*Ta);  // same result
       fprintf(file2,"%f,%12.10f,%12.10f\n",freq,creal(Ta),cimag(Ta));
    }
    fprintf(file2,"END\n");
    fclose(file2);
   return 0;
}

complex double cabl(double freq, double delay, complex double Tin,complex double *ss11, complex double *ss12, double lossf,double dielf)
{ complex double T,Zcab,g,s11,s12,s21,s22,Vin,Iin,Vout,VVin,Z;
  double a,b,d,d2,diel,R,C,L,La,Lb,disp,G;
              { b=0.1175*2.54e-2*0.5; a=0.0362*2.54e-2*0.5; diel = 2.05*dielf; // UT-141C-SP 
                d2 = sqrt(1.0/(PI*4.0*PI*1e-7*5.96e07*0.8*lossf));   // for tinned copper
                d = sqrt(1.0/(PI*4.0*PI*1e-7*5.96e07*lossf));  // skin depth at 1 Hz for copper
              }
  L = (4.0*PI*1e-7/(2.0*PI))*log(b/a);
  C=2.0*PI*8.854e-12*diel/log(b/a);

  La=4.0*PI*1e-7*d/(4.0*PI*a);
  Lb=4.0*PI*1e-7*d2/(4.0*PI*b);
  disp=(La+Lb)/L;
  R = 2.0*PI*L*disp*sqrt(freq*1e6);
  L = L*(1.0+disp/sqrt(freq*1e6));
  G=0;
  if(diel > 1.2) G=2.0*PI*C*freq*1e6*2e-4; // 2e-4 is the loss tangent for teflon
  Zcab = csqrt((I*2*PI*freq*1e6*L+R)/(I*2*PI*freq*1e6*C+G));
  g = csqrt((I*2*PI*freq*1e6*L+R)*(I*2*PI*freq*1e6*C+G));


  T = (50.0-Zcab)/(50.0+Zcab);
  Vin = (cexp(+g*delay*3e08) + T*cexp(-g*delay*3e08));
  Iin = (cexp(+g*delay*3e08) - T*cexp(-g*delay*3e08))/Zcab;
  Vout = (1 + T); // Iout = (1 - T)/Zcab;
  s11 = s22 = ((Vin/Iin) - 50)/((Vin/Iin) + 50);
  VVin = Vin + 50.0*Iin;
  s12 = s21 = (2*Vout/VVin);
  *ss11 = s11;
  *ss12 = s12;

  Z=50.0*(1+Tin)/(1-Tin);
  T = (Z-Zcab)/(Z+Zcab);
  T = T*cexp(-g*2*delay*3e08);
  Z = Zcab*(1+T)/(1-T);
  T = (Z-50.0)/(Z+50.0);
  return T;
}



