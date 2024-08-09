int main(){
    if (fabs(rfi) > 0.0) { // rfi = 2.5 for B18
      file1 = fopen("initial_fourier_model.txt", "w");

      zwt = 1;
      zwtp = 0;
      for (i = 0; i < np; i++)
        if (wtts[i] == 0) wtt[i] = 0;
      for (j = 0; j < 100 && zwt > zwtp; j++) {
        tav = polyfit(pfit, np, data, mcalc, wtt, dataout, fmode);
        rmss = rmscalc(np, dataout, wtt);

        if (j == 0){
          for(i=0;i<np;i++){
            fprintf(file1, "%f %f\n", fstart + i * fstep, dataout[i]);
          }
        }
        float rmsmax = 0;
        float rmsmin = 0;
        float resmax = 0;
        float resmin = 0;
        float zmax = 0;
        float zmin = 0;

        for (i = 0; i < np; i++) {  // sliding rms is better
          a = 0;
          b = 1e-99;
          m = np / 16;
          if (m < 10) m = 10;
          for (k = i - m; k <= i + m; k++) {
            if (k >= 0 && k < np) {
              a += dataout[k] * wtt[k];
              b += wtt[k];
            }
          }
          av = a / b;
          a = 0;
          b = 1e-99;
          for (k = i - m; k <= i + m; k++) {
            if (k >= 0 && k < np) {
              a += (dataout[k] - av) * (dataout[k] - av) * wtt[k];
              b += wtt[k];
            }
          }
          rms = sqrt(a / b);
          b = dataout[i] / (fabs(rfi) * rms);
          if (b > 1) {
            wtt[i] = 0;
            if (nrfi && i + nrfi < np && i - nrfi >= 0)
              for (k = 0; k <= nrfi; k++) wtt[i + k] = wtt[i - k] = 0;
            if (b > 10 && nrfi && i + nrfi * 2 < np && i - nrfi * 2 >= 0)
              for (k = 0; k <= nrfi * 2; k++) wtt[i + k] = wtt[i - k] = 0;
            if (b > 100 && nrfi && i + nrfi * 4 < np && i - nrfi * 4 >= 0)
              for (k = 0; k <= nrfi * 4; k++) wtt[i + k] = wtt[i - k] = 0;
          }
          if(b*rfi < zmin) zmin = b*rfi;
          if(b*rfi > zmax) zmax = b*rfi;
          if(dataout[i] < resmin) resmin = dataout[i];
          if(dataout[i] > resmax) resmax = dataout[i];
          if(rms < rmsmin) rmsmin = rms;
          if(rms > rmsmax) rmsmax = rms;
        }
        zwtp = zwt;
        m = 0;
        for (i = 0; i < np; i++) {
          if (wtt[i] == 0 && fstart + i * fstep > freqstart && fstart + i * fstep < freqstop) m++;
        }
        zwt = m;
        m = 0;
        for (i = 0; i < np; i++) {
          if (fstart + i * fstep > freqstart && fstart + i * fstep < freqstop) m++;
        }
        nnp = m;
        printf("%d rms %f zwt %d res %f %f z %f %f rms %f %f\n", j, rms, zwt, m, np, resmin, resmax, zmin, zmax, rmsmin, rmsmax);
      }
    }
    tav = polyfit(pfit, np, data, mcalc, wtt, dataout, fmode);
}