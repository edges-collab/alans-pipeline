#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void dsmooth(int n, double data[], double data2[], double wt[], double mcalc[], int smooth) {
  double a, b, c;
  int i, j, k;
  for (i = 0; i < n; i++) mcalc[i] = data[i] - data2[i];  // calculate residual
  for (i = 0; i < n; i++) {
    a = b = 0.0;
    for (j = -2.0 * smooth; j <= 2.0 * smooth; j++) {
      k = i + j;
      if (k >= 0 && k < n) {
        c = (k - i) * 2.0 / smooth;
        c = exp(-c * c * 0.69);
        //     if(fabs(k-i)<=smooth/2) c=1; else c=0;  // boxcar
        a += mcalc[k] * wt[k] * c;
        b += c * wt[k];
      }
    }
    if (b > 0.0)
      mcalc[i + n] = a / b;
    else
      mcalc[i + n] = mcalc[i];
  }
  for (i = 0; i < n; i++) {
    //    if(wt[i] > 0.0)
    mcalc[i] = mcalc[i + n] + data2[i];
  }
}

int main(int argc, char *argv[]) {
    FILE *file1;

    file1 = fopen("dsmooth_inputs.txt", "r");
    int n = 0;
    int smooth = 8;
    double tmp;
    double data[3000];
    double data2[3000];
    double wt[3000];
    double mcalc[6000];

    while (fscanf(file1, "%le %le %le", &data[n], &data2[n], &wt[n]) != EOF) {
        n++;
    }
    fclose(file1);
    
    
    dsmooth(n, data, data2, wt, mcalc, smooth);

    file1 = fopen("dsmooth_outputs.txt", "w");
    for (int i = 0; i < n; i++) {
        fprintf(file1, "%1.10e %1.10e\n", mcalc[i] - data2[i], mcalc[i]);
    }

    printf("Done, and %d\n", 2 & 1);
}