#include <complex.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/io.h>
#include <time.h>
#define PI 3.1415926536
#define TWOPI 6.28318530717958
#define NDATA 8192  // was 32768
#define NSPEC 2000  // was 500

void outspec(int *, double *, double *, double *, double *, double *, int, int, char *);
void plotfspec(int *, double *, double *, double *, int, double, char *, char *, int, double, int, char *);
void plotfspec2(int *, double *, double *, double *, int, double, double);
void zplott(int, double *, double *, double *, double *, double *, int, char *);
void zplott2(int, double *, double *, double *, double *, double *, double *, int, char *, char *);
double fition(int, int, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, double,
              int, int, double, double);
double polyfitr(int, int, double *, double *, double *, double *, double *, long double *, double *);
double fitfun(int, int, int, double *);
void qrd(long double *, int, double *);
double tosecs(int, int, int, int, int);
void radec_azel(double, double, double, double *, double *);
void moonradec(double, double *, double *, double, double);
void sunradec(double, double *, double *);
int moon(double, double);
int sun(double, double);
double gst(double);
double rmscalc(int, double *, double *, double *);
void eorsnr(int, double *, double *, double *, double *, double, double, double, int, int, double, double, int, double, double, char *);
double gauss();
double agauss(double, double, double);
void fit(int, int, int, double *, double *, double *, double *, double *, double *, double, double, double);
double chifit(int, int, double *, double *, double *, double *, int, double, double *, double *, double *, double *, double *, int, double, double,
              unsigned char *, double, double, double, double);
int chisolve(int, double *, double *, double *, double *, int, int, double, double, double, double);
void mapp(unsigned char *, double, double, double, double, double, double, double, double, double, double, double, double, double);
void plotq(unsigned char *, double, double, double, double, double, double, double, double, char *);

static double fitf[NDATA * 40], bbrr[NDATA], mcalc[NDATA * 40];
static long double aarr[NDATA];

int main(int argc, char *argv[]) {
  char name[255], buf[32768], ti[255], tti[255], txt[255];
  static char title[100 * NSPEC], title2[100 * NSPEC], title3[100 * NSPEC], title4[100 * NSPEC], titl[100 * NSPEC];
  static char info[100 * 1000], info2[100 * NSPEC], info3[100 * NSPEC], info4[256];
  FILE *file3;
  static double spec[NSPEC * NDATA], spdiff[NSPEC * NDATA], freq[NSPEC * NDATA], wtt[NSPEC * NDATA], specout[NSPEC * NDATA], wtav[NDATA],
      freqpl[NSPEC * NDATA];
  static double sppl[NSPEC * NDATA], spplm[NSPEC * NDATA], wtpl[NSPEC * NDATA], t150p[NSPEC];
  static double dnav[NDATA], dnavm[NDATA], wtdn[NDATA], timm[NDATA];
  static int typ[NDATA], typp[NDATA];
  double dmax, frq, spe, frqp, mde, wt;
  double rms, maxrms, t150, f150, spind, spcurv, ionabs, ionemm, speor, cov, lim, limm, fstart, fstop, wid, eora, aeora, seor, feor, rfi, mchk, tchk,
      schk, test;
  double snr, adnoise, gha, t150dn, diff, avrms, tau, fsmax, fsmin;
  int i, j, m, yr, day, hr, min, sc, nplot, j3, nfit, nn[NSPEC], nnpl[NSPEC], out, wtmode, sig, ndn, g10, alt, half, sim, imode, tmode, pmode, rr,
      resid, sm, md;
  dmax = 10;
  lim = 1e6;
  fstart = 0;
  fstop = 200;
  nfit = 4;
  seor = 0;
  adnoise = 0;
  out = 0;
  sig = 0;
  g10 = 1;
  diff = alt = half = rfi = mchk = schk = imode = tmode = pmode = rr = test = resid = 0;
  eora = aeora = 0;
  feor = 150;
  wid = 10;
  wtmode = 0;
  sim = 0;
  tchk = -1e99;
  tau = 0;
  sm = 0;
  md = 0;
  fsmax = 200;
  fsmin = 40;
  sprintf(ti, "GHA");
  sscanf(argv[1], "%s", name);
  for (m = 0; m < argc; m++) {
    sscanf(argv[m], "%79s", buf);
    if (strstr(buf, "-dmax")) { sscanf(argv[m + 1], "%lf", &dmax); }
    if (strstr(buf, "-lim")) { sscanf(argv[m + 1], "%lf", &lim); }
    if (strstr(buf, "-fstart")) { sscanf(argv[m + 1], "%lf", &fstart); }
    if (strstr(buf, "-fstop")) { sscanf(argv[m + 1], "%lf", &fstop); }
    if (strstr(buf, "-fsmax")) { sscanf(argv[m + 1], "%lf", &fsmax); }
    if (strstr(buf, "-fsmin")) { sscanf(argv[m + 1], "%lf", &fsmin); }
    if (strstr(buf, "-eor")) { sscanf(argv[m + 1], "%lf", &eora); }
    if (strstr(buf, "-aeor")) { sscanf(argv[m + 1], "%lf", &aeora); }
    if (strstr(buf, "-feor")) { sscanf(argv[m + 1], "%lf", &feor); }
    if (strstr(buf, "-seor")) { sscanf(argv[m + 1], "%lf", &seor); }
    if (strstr(buf, "-wid")) { sscanf(argv[m + 1], "%lf", &wid); }
    if (strstr(buf, "-nfit")) { sscanf(argv[m + 1], "%d", &nfit); }
    if (strstr(buf, "-wtmode")) { sscanf(argv[m + 1], "%d", &wtmode); }
    if (strstr(buf, "-imode")) { sscanf(argv[m + 1], "%d", &imode); }
    if (strstr(buf, "-pmode")) { sscanf(argv[m + 1], "%d", &pmode); }
    if (strstr(buf, "-out")) { sscanf(argv[m + 1], "%d", &out); }
    if (strstr(buf, "-sig")) { sscanf(argv[m + 1], "%d", &sig); }
    if (strstr(buf, "-test")) { sscanf(argv[m + 1], "%lf", &test); }
    if (strstr(buf, "-tau")) { sscanf(argv[m + 1], "%lf", &tau); }
    if (strstr(buf, "-adnoise")) { sscanf(argv[m + 1], "%lf", &adnoise); }
    if (strstr(buf, "-rfi")) { sscanf(argv[m + 1], "%lf", &rfi); }
    if (strstr(buf, "-mchk")) { sscanf(argv[m + 1], "%lf", &mchk); }
    if (strstr(buf, "-tchk")) { sscanf(argv[m + 1], "%lf", &tchk); }
    if (strstr(buf, "-schk")) { sscanf(argv[m + 1], "%lf", &schk); }
    if (strstr(buf, "-g10")) g10 = 10;
    if (strstr(buf, "-gg100")) g10 = 100;
    if (strstr(buf, "-diff")) { sscanf(argv[m + 1], "%lf", &diff); }
    if (strstr(buf, "-rr")) { sscanf(argv[m + 1], "%d", &rr); }
    if (strstr(buf, "-alt")) alt = 1;
    if (strstr(buf, "-half")) half = 1;
    if (strstr(buf, "-sim")) { sscanf(argv[m + 1], "%d", &sim); }
    if (strstr(buf, "-resid")) { sscanf(argv[m + 1], "%d", &resid); }
    if (strstr(buf, "-sm")) { sscanf(argv[m + 1], "%d", &sm); }
    if (strstr(buf, "-md")) { sscanf(argv[m + 1], "%d", &md); }  // plot signature
    if (strstr(buf, "-ti")) {
      sscanf(argv[m + 1], "%s", ti);
      tmode = 1;
    }
    if (strstr(buf, "-txt")) sscanf(argv[m + 1], "%s", txt);
    if (strstr(buf, "-date")) { tmode = 2; }
    if (strstr(buf, "-datt")) { tmode = 3; }
  }
  if ((file3 = fopen(name, "r")) == NULL) {
    printf("%s error\n", name);
    return 0;
  }
  if (rr) srand(rr);
  i = 0;
  frqp = 1e99;
  j = -1;
  if (fstop < 150 && feor > 100) feor = 75;
  while (fgets(buf, 32768, file3) != 0) {
    if (strstr(buf, ti)) sscanf(buf, "%*s %s", tti);
    if (strstr(buf, "freq") && !strstr(buf, "nan")) {
      // freq 100.0045 tantenna 746.8551 K skymodel 746.9250 K wt 1 2015:227:20:43:03
      sscanf(buf, "%*s %lf %*s %lf %*s %*s %lf %*s %*s %lf %d:%d:%d:%d:%d", &frq, &spe, &mde, &wt, &yr, &day, &hr, &min, &sc);
      if (frq < frqp) {
        i = 0;
        j++;  // printf("j %d freq %f %f\n",j,frq,frqp);
        sprintf(&titl[j * 100], "%2d:%03d", yr - 2010, day);
        if (tmode == 3) sprintf(&titl[j * 100], "%2d:%03d:%02d", yr - 2010, day, hr);
        info2[(j + 1) * 100] = 0;
        sscanf(buf, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %s %s", &title[j * 100], &title3[j * 100]);
        if (tmode == 1) sprintf(&title3[j * 100], "%s%s", ti, tti);
      }
      frqp = frq;
      freq[i + j * NDATA] = frq;
      if (frq < fstart || frq > fstop) wt = 0;
      if (sim) { // not done for B18
        spe = sim * 300.0 * pow(frq / 150, -2.6 + 0.1 * log(frq / 50) + 0.0 * pow(log(frq / 75), 2.0)) * exp(-2e-2 * 150.0 * 150.0 / (frq * frq)) +
              2.0e-2 * 3e3 * pow(frq / 150, -2.0);
        //          if(sim) {spe = sim *
        //          300.0*pow(frq/150,-2.6+0.1*log(frq/50)+0.0*pow(log(frq/75),2.0))*exp(-2e-3*150.0*150.0/(frq*frq))+2.0e-3*3e3*pow(frq/150,-2.0);
        //                   spe += sim *
        //                   300.0*pow(frq/150,-2.6+0.3*log(frq/50)+0.1*pow(log(frq/75),2.0))*exp(-2e-2*150.0*150.0/(frq*frq))+2.0e-2*3e3*pow(frq/150,-2.0);
      }
      // Neither of the following done for B18
      if (eora && (sig % 10) == 0) spe += -eora * exp(-0.693 * (frq - feor) * (frq - feor) / (wid * wid * 0.25));  // negative for absorption
      if (aeora && (sig % 10) == 0) spe += -aeora * agauss(frq, feor, wid);                                        // negative for absorption
      //          if(test) spe += -test*(1-exp(-flt*exp(-0.693*(frq-feor)*(frq-feor)/(wid*wid*0.25))));  // test positive for absorption
      if (test) // not done for B18
        spe += -test * (1 - exp(-tau * exp(-(frq - feor) * (frq - feor) * (-log(-log((1 + exp(-tau)) / 2) / tau)) / (wid * wid * 0.25)))) /
               (1 - exp(-tau));  // true width
      if (eora && (sig % 10) == 1) spe += eora * 0.5 * (tanh((1420 / frq - 1420 / feor) / (wid * 1420 / (feor * feor))) + 1.0);
      if (adnoise) {
        spe += adnoise * 1e-3 * gauss();
        if (frq >= fstart && frq <= fstop) wt = 1;
      }
      if (resid) // not done for B18
        spec[i + j * NDATA] = spe - mde;
      else
        spec[i + j * NDATA] = spe;
      wtt[i + j * NDATA] = wt;
      timm[j] = tosecs(yr, day, hr, min, sc);
      typ[j] = 1;
      if (alt) typ[j] = j % 2;
      i++;
    }
    nn[j] = i;
    if (strstr(buf, ti)) sscanf(buf, "%*s %s", &info2[(j + 1) * 100]);
  }
  fclose(file3);
  nplot = j + 1;
  if (half && !alt) {
    for (j = 0; j < nplot; j++)
      if (j < nplot / 2)
        typ[j] = 0;
      else
        typ[j] = 1;
  }
  if (half && alt) {
    for (j = 0; j < nplot; j++)
      if (j < nplot / 2)
        typ[j] = 1;
      else
        typ[j] = 0;
  }
  maxrms = -1e99;
  j3 = 0;
  avrms = 0;
  for (j = 0; j < nplot; j++) {
    for (i = 0; i < nn[j]; i++) {
      spdiff[i] = spec[i + j * NDATA];
      if (j < nplot - 1)
        m = j + 1;
      else
        m = j;
      // Not done in B18:
      if (diff < 0.0) spdiff[i] = spec[i + m * NDATA] - spec[i + j * NDATA];  // straight difference
      wtav[i] = wtt[i + j * NDATA];
    }
    f150 = fition(nn[j], abs(nfit), freq, spdiff, wtav, specout, &t150, &spind, &spcurv, &ionabs, &ionemm, &speor, &cov, feor, wid, sig, 0, 0, 0);
    rms = rmscalc(nn[j], spdiff, specout, wtav);
    for (i = 0; i < nn[0]; i++)
      if (rfi && (spdiff[i] - specout[i]) > rms * rfi) wtav[i] = 0;
    if (rms > maxrms) maxrms = rms;
    limm = lim;
    if (typ[j] == 0) limm = 0;
    m = 0;
    printf(">> tchk %f t150 %f sig %d\n", tchk, t150, sig);

    // For B18: mchk = 0, schk = 0 (so they both pass). t150 is set in the fit above
    // and is, I think actually T75 modelled.
    if (rms < limm && moon(timm[j], mchk) && sun(timm[j], schk) && t150 > tchk) {
      for (i = 0; i < nn[j]; i++) {
        sppl[i + j3 * NDATA] = spdiff[i] - specout[i];
        spplm[i + j3 * NDATA] = specout[i];
        // sppl is a residual - tests show it can be a total and it doesn't effect difference
        wtpl[i + j3 * NDATA] = wtav[i];
        freqpl[i + j3 * NDATA] = freq[i];
        nnpl[j3] = nn[j];
        t150p[j3] = t150;
      }
      sprintf(&title2[j3 * 100], "%s %d", &titl[j * 100], typ[j]);
      if (tmode == 3) sprintf(&title2[j3 * 100], "%s ", &titl[j * 100]);
      sprintf(&title4[j3 * 100], "%s", &title[j * 100]);
      if (info2[j * 100]) sprintf(&title2[j3 * 100], "%s %s %d", ti, &info2[j * 100], typ[j]);
      if (title3[j * 100] && tmode != 2 && tmode != 3)
        sprintf(&title2[j3 * 100], " %s", &title3[j * 100]);  // needs condx to avoid getting nline instead of date
      sprintf(&info[j3 * 100], "rms %3.1e", rms);
      avrms += rms;
      typp[j3] = typ[j];
      m = 1;
      j3++;
    }
    printf("t%d %8.3f rms %8.5f %s %d spind %10.3f ionabspercent %7.2f ionemmK %7.0f spcurv %10.3f cov %10.3f m %d nline %s\n", (int)f150, t150, rms,
           &title[j * 100], typ[j], spind, ionabs, ionemm, spcurv, cov, m, &title3[j * 100]);
  }
  if (!j3) {
    printf("no data accepted\n");
    return 0;
  }
  for (i = 0; i < nn[0]; i++) dnav[i] = dnavm[i] = wtdn[i] = 0;
  ndn = t150dn = 0;
  for (j = 0; j < j3; j++) {
    if (typp[j] == 1) {
      t150dn += t150p[j];
      ndn++;
    }
    for (i = 0; i < nn[j]; i++) {
      if (typp[j] == 1) {
        dnav[i] += wtpl[i + j * NDATA] * sppl[i + j * NDATA];
        wtdn[i] += wtpl[i + j * NDATA];
        dnavm[i] += spplm[i + j * NDATA];
      }
    }
  }
  for (i = 0; i < nn[0]; i++) {
    if (wtdn[i]) dnav[i] = dnav[i] / wtdn[i];
  }
  //      for(i=0;i<nn[0];i++) spdiff[i]= dnav[i];
  for (i = 0; i < nn[0]; i++) spdiff[i] = dnav[i] + dnavm[i] / ndn;
  for (i = 0; i < nn[0]; i++) {
    if (wtdn[i] < ndn / 3)
      wtav[i] = 0;
    else
      wtav[i] = 1;
  }
  //      for(i=0;i<nn[0];i++) wtav[i]=1;
  for (i = 0; i < nn[0]; i++) {
    if (wtdn[i] == 0) wtav[i] = 0;
  }
  if (diff > 0)
    for (i = 0; i < nn[0]; i++) {
      spdiff[i] = (sppl[i + NDATA] - diff * sppl[i]) / (1.0 - diff);
      // note correction for amplitude as per memo 247
      if (diff == 1) spdiff[i] = (sppl[i + NDATA] - sppl[i]);
      //      if(diff > 0) for(i=0;i<nn[0];i++) {spdiff[i] = (sppl[i+NDATA] + spplm[i+NDATA] - diff*(sppl[i]+spplm[i]))/(1.0-diff); // tested made no
      //      diff if(diff > 0) for(i=0;i<nn[0];i++) {spdiff[i] = 4000.0*(sppl[i+NDATA] + spplm[i+NDATA])/(sppl[i]+spplm[i]); // ratio test
      //            if(i==nn[0]/2) printf("rmsin sp %f %f\n",spdiff[i],spplm[i+NDATA]);
      if (wtpl[i + NDATA] && wtpl[i])
        wtdn[i] = 1;
      else
        wtdn[i] = 1;
    }
  if (rfi) {
    fition(nn[0], 7, freq, spdiff, wtav, specout, &t150, &spind, &spcurv, &ionabs, &ionemm, &speor, &cov, feor, wid, 10, 0, 0,
           tau);  // rfi filt with 7 term poly
    rms = rmscalc(nn[0], spdiff, specout, wtav);
    for (i = 0; i < nn[0]; i++)
      if ((spdiff[i] - specout[i]) > rms * rfi) wtav[i] = 0;
  }
  if (nfit > 0)
    fition(nn[0], nfit, freq, spdiff, wtav, specout, &t150, &spind, &spcurv, &ionabs, &ionemm, &speor, &cov, feor, wid, sig, imode,
           t150dn / (ndn + 1e-6), 0);
  else
    fition(nn[0], -abs(fabs(nfit) + 1), freq, spdiff, wtav, specout, &t150, &spind, &spcurv, &ionabs, &ionemm, &speor, &cov, feor, wid, sig, 0, 0,
           tau);
  rms = rmscalc(nn[0], spdiff, specout, wtav);
  sprintf(&title2[j3 * 100], "%s ", " av");
  if (diff) sprintf(&title2[j3 * 100], "%s ", " diff");
  sprintf(&info[j3 * 100], "rms %3.1e scale x %d", rms, g10);
  sprintf(info3, "%s %s%02.0f", titl, ti, gha);
  sprintf(info4, "%s", ti);
  //        sprintf(&info[j3*100],"rms %4.0f %4.0f/div mK",rms*1e3,2e3*dmax/sqrt(j3));
  if (sm) {
    for (i = 0; i < nn[0]; i++) spdiff[i] = spdiff[i] - specout[i];
    fition(nn[0], sm, freq, spdiff, wtav, specout, &t150, &spind, &spcurv, &ionabs, &ionemm, &speor, &cov, feor, wid, sig, imode,
           t150dn / (ndn + 1e-6), 0);
    for (i = 0; i < nn[0]; i++) {
      spdiff[i] = specout[i];
      specout[i] = 0;
    }
    rms = rmscalc(nn[0], spdiff, specout, wtav);
    sprintf(&info[j3 * 100], "rms %3.1e scale x %d", rms, g10);
  }
  for (i = 0; i < nn[0]; i++) {
    sppl[i + j3 * NDATA] = spdiff[i] - specout[i];
    wtpl[i + j3 * NDATA] = wtav[i];
    freqpl[i + j3 * NDATA] = freq[i];
  }
  //      if(imode) for(i=0;i<nn[0];i++) sppl[i+j3*NDATA]=spdiff[i]-specout[i];  // for plot
  nnpl[j3] = nn[0];
  j3++;
  snr = speor / (sqrt(fabs(cov)) * rms);
  printf("nplot %d j3 %d ndn %d\n", nplot, j3, ndn);
  printf(
      "t%d %7.3f maxrms %5.1f avrms %6.4f rmsofav %6.4f K j3 %d %s t150dn %6.4f spind %6.4f ionabs %6.4f ionemm %6.4f spcurv %6.4f eor_K %6.4f snr "
      "%f\n",
      (int)feor, t150, maxrms, avrms / (j3 - 1.0), rms, j3, &title4[(j3 - 2) * 100], t150dn / (ndn + 1e-6), spind, ionabs, ionemm, spcurv, speor,
      snr);
  if (j3 > 2) plotfspec(nnpl, freqpl, sppl, wtpl, j3, dmax, title2, info, g10, avrms / (j3 - 1.0), pmode, txt);
  if (j3 == 2) plotfspec2(nnpl, freqpl, sppl, wtpl, j3, dmax, rms);
  if (seor < 0) eorsnr(nn[0], freq, spdiff, wtav, specout, fstart, fstop, wid, sig, nfit, seor, tau, md, fsmax, fsmin, info4);
  if (out) outspec(nnpl, freqpl, sppl, spdiff, specout, wtpl, j3, out, info3);
  if (seor == -3) chisolve(nn[0], freq, spdiff, wtav, specout, sig, nfit, fstart, fstop, fsmax, fsmin);
  return 0;
}

void plotfspec(int nn[], double freqq[], double data[], double wt[], int nplot, double ddmax, char title[], char info[], int g10, double avrms,
               int pmode, char txxt[])
// plot the spectrum
{
  char txt[256];
  int k, iter, np, nnplot;
  double x, y, z, xp, yp, ypp, dmax, dmin, f, totpp, scale, step, h, s, b;
  double xoffset, yoffset, fstart, fstop, dd;
  FILE *file;

  sprintf(txt, "spe.pos");
  if ((file = fopen(txt, "w")) == NULL) {
    printf("cannot open %s:\n", txt);
    return;
  }
  fprintf(file, "%%!PS-Adobe-3.0 EPSF-3.0\n%c%cBoundingBox:  0 0 612 700\n%c%cEndProlog\n", '%', '%', '%', '%');
  fprintf(file, "1 setlinewidth\n");
  fprintf(file, "/Times-Roman findfont\n 14 scalefont\n setfont\n");

  fstart = fstop = 0;
  xoffset = 50.0;
  yoffset = 100.0;
  dmin = -ddmax;
  dmax = ddmax;
  scale = dmax - dmin;
  for (y = 0; y < 2; y++)
    fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset, y * 480 + yoffset, xoffset + 400.0,
            y * 480 + yoffset);
  fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset, yoffset, xoffset, 480.0 + yoffset);
  fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset + 400.0, yoffset, xoffset + 400.0,
          480.0 + yoffset);
  fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset + 0.0, yoffset, xoffset + 0.0,
          480.0 + yoffset);
  if (pmode == 1) {
    nnplot = nplot - 1;
    dd = 480.0;
  } else {
    nnplot = nplot;
    dd = 440.0;
  }
  for (iter = 0; iter < nnplot; iter++) {
    np = nn[iter];
    fstart = freqq[0];
    fstop = freqq[np - 1 + iter * NDATA];
    fstop = (int)(fstop + 0.5);  // set to nearest MHz
    yp = 0;
    xp = 0;
    totpp = 0;
    ypp = 0;
    for (k = 0; k < np; k++) {
      h = s = b = 0;
      x = (freqq[k + iter * NDATA] - fstart) * 400.0 / (fstop - fstart);
      totpp = (data[k + iter * NDATA] - dmin) / scale;
      y = 480.0 - dd / nnplot - iter * dd / nnplot + totpp * dd / nnplot;
      z = (480.0 - 440.0 / nplot - (nplot - 2) * 440.0 / nplot) / 2.0;
      if (z < 20) z = 20;
      if (iter == nplot - 1) y = z + g10 * (data[k + iter * NDATA] / scale) * 440.0 / nplot;
      if (y > 480.0) y = 480.0;
      if (y < 0.0) y = 0.0;
      if (k == 0 || x < xp) {
        xp = x;
        yp = y;
      }
      h = s = b = 0;
      if (k && wt[k + iter * NDATA] > 0.0 && wt[k - 1 + iter * NDATA] > 0.0) {
        fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n %5.2f %5.2f %5.2f sethsbcolor stroke\n", xp + xoffset, yp + yoffset,
                x + xoffset, y + yoffset, h, s, b);
        ypp = yp;
      }
      xp = x;
      yp = y;
    }
    totpp = -dmin / scale;
    if (iter < nplot - 1)
      ypp = 480.0 - dd / nnplot - iter * dd / nnplot + totpp * dd / nnplot;
    else
      ypp = 20;
    sprintf(txt, "%s %s", &title[iter * 100], &info[iter * 100]);
    //        if(nplot < 25 || iter==nplot-1) fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 400.0,ypp + yoffset,txt);
    if (nplot > 50 && iter < nplot - 1) fprintf(file, "/Times-Roman findfont\n 7 scalefont\n setfont\n");
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 400.0, ypp + yoffset, txt);
    fprintf(file, "/Times-Roman findfont\n 14 scalefont\n setfont\n");
  }

  step = 5.0;
  if (fstop - fstart > 60) step = 10;
  if (fstop - fstart > 100) step = 20;
  for (f = fstart; f <= fstop + step * 0.1; f += step) {
    x = (f - fstart) * 400.0 / (fstop - fstart);
    y = 0;
    fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, y + yoffset, x + xoffset,
            y + 5.0 + yoffset);
    //        sprintf(txt, "%5.1f", f);
    sprintf(txt, "%4.0f", f);
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 15.0, y - 15.0 + yoffset, txt);
  }

  for (iter = 0; iter < nnplot; iter++) {
    for (f = dmin; f <= dmax; f += (dmax - dmin)) {
      x = 0;
      y = (f - dmin) / (dmax - dmin);
      if (iter < nplot - 1) {
        y = 480.0 - dd / nnplot - iter * dd / nnplot + y * dd / nnplot;
        //        else y = 20;
        fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, y + yoffset, x + xoffset + 5,
                y + yoffset);
      }
      //           {sprintf(txt, "%4.0f K", f); fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n",x + xoffset - 45.0,y - 1.0+yoffset, txt);}
    }
  }
  if (pmode == 2) {
    fprintf(file, "0.1 setlinewidth\n");
    step = 5.0;
    z = (480.0 - 440.0 / nplot - (nplot - 2) * 440.0 / nplot) / 2.0;
    if (fstop - fstart > 60) step = 10;
    if (fstop - fstart > 100) step = 20;
    for (f = fstart; f <= fstop + step * 0.1; f += step) {
      x = (f - fstart) * 400.0 / (fstop - fstart);
      y = 0;
      fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, y + z - 0.5 * dd / nplot + yoffset,
              x + xoffset, y + z + 0.5 * dd / nplot + yoffset);
      sprintf(txt, "%4.0f", f);
      fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 15.0, y - 15.0 + yoffset, txt);
    }

    sprintf(txt, "%5.2f", scale * 0.5);
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset - 35.0, yoffset - 4.0 + z + 0.5 * dd / nplot, txt);

    sprintf(txt, "%5.2f", 0.0);
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset - 35.0, yoffset - 4.0 + z + 0.0 * dd / nplot, txt);

    sprintf(txt, "%5.2f", -scale * 0.5);
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset - 35.0, yoffset - 4.0 + z - 0.5 * dd / nplot, txt);

    for (iter = 0; iter < 5; iter++) {
      x = 0;
      y = z + (iter - 2) * 0.25 * dd / nplot;
      fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, y + yoffset, x + xoffset + 400.0,
              y + yoffset);
    }
    fprintf(file, "1 setlinewidth\n");
  }
  /*
      for (f = dmin; f <=dmax; f += (dmax - dmin)/4.0) {
          x = 0;
          z = (480.0 - 440.0/nplot - (nplot-2)*440.0/nplot)/2.0;
          if(z < 20) z = 20;
          y = z + 4.0*(f/scale)*440.0/nplot;
          fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n",
             x + xoffset,y+yoffset,x + xoffset + 10,y+yoffset);
      }
  */
  sprintf(txt, " temperature %6.2f K  per_division", dmax - dmin);
  //    fprintf(file, "%6.2f %6.2f moveto\n 90 rotate\n (%s) show\n -90 rotate\n",xoffset -52.0, 265.0, txt);
  fprintf(file, "%6.2f %6.2f moveto\n 90 rotate\n (%s) show\n -90 rotate\n", xoffset - 20.0, 230.0, txt);

  sprintf(txt, "Frequency (MHz)");
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 160.0, 70.0, txt);

  sprintf(txt, "avrms %6.4f", avrms);
  if (pmode == 0) fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 460.0, 50.0, txt);
  if (txxt[0]) fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 400.0, 50.0, txxt);

  fprintf(file, "showpage\n%c%cTrailer\n", '%', '%');
  fclose(file);
}

void plotfspec2(int nn[], double freqq[], double data[], double wt[], int nplot, double ddmax, double rms)
// plot the spectrum
{
  char txt[256];
  int k, iter, np;
  double x, y, xp, yp, dmax, dmin, f, totpp, scale, step, h, s, b;
  double xoffset, yoffset, fstart, fstop;
  FILE *file;

  sprintf(txt, "spe.pos");
  if ((file = fopen(txt, "w")) == NULL) {
    printf("cannot open %s:\n", txt);
    return;
  }
  fprintf(file, "%%!PS-Adobe-3.0 EPSF-3.0\n%c%cBoundingBox:  0 0 612 700\n%c%cEndProlog\n", '%', '%', '%', '%');
  fprintf(file, "1 setlinewidth\n");
  fprintf(file, "/Times-Roman findfont\n 14 scalefont\n setfont\n");

  fstart = fstop = 0;
  xoffset = 50.0;
  yoffset = 100.0;
  dmin = -ddmax;
  dmax = ddmax;
  scale = dmax - dmin;
  for (y = 0; y < 2; y++)
    fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset, y * 480 + yoffset, xoffset + 400.0,
            y * 480 + yoffset);
  fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset, yoffset, xoffset, 480.0 + yoffset);
  fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset + 400.0, yoffset, xoffset + 400.0,
          480.0 + yoffset);
  fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset + 0.0, yoffset, xoffset + 0.0,
          480.0 + yoffset);
  for (iter = 1; iter < nplot; iter++) {
    np = nn[iter];
    fstart = freqq[0];
    fstop = freqq[np - 1 + iter * NDATA];
    fstop = (int)(fstop + 0.5);  // set to nearest MHz
    yp = 0;
    xp = 0;
    totpp = 0;

    for (k = 0; k < np; k++) {
      h = s = b = 0;
      x = (freqq[k] - fstart) * 400.0 / (fstop - fstart);
      totpp = (data[k + NDATA] - dmin) / scale;
      y = totpp * 480.0;
      if (y > 480.0) y = 480.0;
      if (y < 0.0) y = 0.0;
      if (k == 0) {
        xp = x;
        yp = y;
      }
      h = s = b = 0;
      if (k && wt[k + iter * NDATA] > 0.0 && wt[k - 1 + iter * NDATA] > 0.0) {
        fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n %5.2f %5.2f %5.2f sethsbcolor stroke\n", xp + xoffset, yp + yoffset,
                x + xoffset, y + yoffset, h, s, b);
      }
      xp = x;
      yp = y;
    }
  }
  if (fstop < 100)
    step = 5.0;
  else
    step = 10.0;
  fprintf(file, "0.1 setlinewidth\n");
  for (f = fstart; f <= fstop + step * 0.1; f += step) {
    x = (f - fstart) * 400.0 / (fstop - fstart);
    y = 0;
    fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, y + yoffset, x + xoffset,
            y + 480.0 + yoffset);
    sprintf(txt, "%4.0f", f);
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 15.0, y - 15.0 + yoffset, txt);
  }

  for (f = dmin; f <= dmax; f += (dmax - dmin) * 0.05) {
    x = 0;
    y = (f - dmin) / (dmax - dmin);
    y = 480.0 * y;
    fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, y + yoffset, x + xoffset + 400.0,
            y + yoffset);
    if (dmax < 10)
      sprintf(txt, "%3.1f", f);
    else
      sprintf(txt, "%3.0f", f);
    if (dmax < 1) sprintf(txt, "%3.0f", f * 1e3);
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 26.0, y - 0.0 + yoffset, txt);
  }
  fprintf(file, "1 setlinewidth\n");
  sprintf(txt, " temperature (K)");
  if (dmax < 1) sprintf(txt, " temperature (mK)");
  fprintf(file, "%6.2f %6.2f moveto\n 90 rotate\n (%s) show\n -90 rotate\n", xoffset - 30.0, 305.0, txt);
  sprintf(txt, "rms %3.0f mK", rms * 1e3);
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 410.0, 150.0, txt);
  sprintf(txt, "Frequency (MHz)");
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 160.0, 70.0, txt);
  fprintf(file, "showpage\n%c%cTrailer\n", '%', '%');
  fclose(file);
}

/* Convert to Seconds since New Year 1970 */
double tosecs(int yr, int day, int hr, int min, int sec) {
  int i;
  double secs;
  secs = (yr - 1970) * 31536000.0 + (day - 1) * 86400.0 + hr * 3600.0 + min * 60.0 + sec;
  for (i = 1970; i < yr; i++) {
    if ((i % 4 == 0 && i % 100 != 0) || i % 400 == 0) secs += 86400.0;
  }
  if (secs < 0.0) secs = 0.0;
  return secs;
}

/* Calculate Sun ra and dec (approximate) */
/* see Astronomical Almanac page C24 Sun 1999 */
void sunradec(double time, double *ra, double *dec) {
  double n, g, lonn, ecl;
  n = -365.5 + (time - tosecs(1999, 1, 0, 0, 0)) / 86400.0;
  g = (357.528 + 0.9856003 * n) * PI / 180.0;
  lonn = (280.460 + 0.9856474 * n + 1.915 * sin(g) + 0.02 * sin(2 * g)) * PI / 180.0;
  ecl = (23.439 - 0.0000004 * n) * PI / 180.0;
  *ra = atan2(sin(lonn) * cos(ecl), cos(lonn));
  if (*ra < 0) *ra += TWOPI;
  *dec = asin(sin(lonn) * sin(ecl));
}

void moonradec(double time, double *ra, double *dec, double lat, double lon) {
  double ttime, asnode, amon, peri, em, aim, aam, vsm, alamm, moonsun, evn, var;
  double x, y, z, xx, yy, zz, ram, decm, ha, inc;
  /* calc moon ra and dec  */
  /* see notes and formulae Astronomical Almanac page D2 Moon, 1999 */
  ttime = tosecs(1999, 0, 0, 0, 0); /* Jan 0 1999 */
  ttime = (time - ttime) / 86400.00;
  /* asnode=long of mean ascending node */
  asnode = 144.452077 - 0.05295377 * ttime;
  /* amon=omga plus mean lunar longitude */
  amon = 69.167124 + 13.17639648 * ttime;
  /* peri=asnode plus mean lunar longitude of perigee */
  peri = 42.524057 + 0.11140353 * ttime;
  /* moonsun is the elongation of moon from the sun */
  moonsun = 149.940812 + 12.19074912 * ttime;
  /* em is the eccentricity of lunar orbit */
  em = 0.054900489;
  /* aim=inclination of lunar orbit to ecliptic */
  aim = 5.1453964 * PI / 180.0;
  /* vsm=true anomaly */
  /* the following are correction terms */
  vsm = 2.0 * em * sin((amon - peri) * PI / 180.0);                       /* elliptical orbit */
  evn = (1.274 / 57.3) * sin((2 * moonsun - (amon - peri)) * PI / 180.0); /* evection */
  var = (0.658 / 57.3) * sin(2 * moonsun * PI / 180.0);                   /* variation */
  alamm = (amon - asnode) * PI / 180.0 + vsm + evn + var;
  x = cos(alamm);
  y = sin(alamm);
  z = 0;
  xx = x;
  yy = y * cos(aim);
  zz = y * sin(aim);
  ram = atan2(yy, xx) + asnode * PI / 180.0;
  decm = atan2(zz, sqrt(xx * xx + yy * yy));
  x = cos(ram) * cos(decm);
  y = sin(ram) * cos(decm);
  z = sin(decm);
  inc = 23.45 * PI / 180.0;
  xx = x;
  yy = y * cos(inc) - z * sin(inc);
  zz = z * cos(inc) + y * sin(inc);
  /* aam is the semi-major axis of orbit earth radii */
  aam = 60.2665;
  z = zz - sin(lat) / aam; /* correct for parallax */
  ha = gst(time) + lon;    // longitude east
  x = xx - cos(lat) * cos(ha) / aam;
  y = yy - cos(lat) * sin(ha) / aam;
  *ra = atan2(y, x);
  if (*ra < 0) *ra += TWOPI;
  *dec = atan2(z, sqrt(x * x + y * y));
}

double gst(double ttime) {
  double secs, pdum;
  int i;
  secs = (2011 - 1970) * 31536000.0 + 17.0 * 3600.0 + 15.0 * 60.0 + 58.0778;
  for (i = 1970; i < 2011; i++) {
    if ((i % 4 == 0 && i % 100 != 0) || i % 400 == 0) secs += 86400.0;
  }

  return (modf((ttime - secs) / 86164.09053, &pdum) * TWOPI);
  // 17 15 58.0778 UT at 0hr newyear2011
}

void radec_azel(double ha, double dec, double latt, double *azs, double *elevs) {
  /* convert from sky to antenna coords (azel mount) */
  /* input: ha,dec,latt
   *    output: azs=azimuth of source
   *       elevs=elevation of source
   *        */
  double p, w, r, zen, north;
  p = sin(dec);
  w = sin(ha) * cos(dec);
  r = cos(ha) * cos(dec);
  zen = r * cos(latt) + p * sin(latt);
  north = -r * sin(latt) + p * cos(latt);
  *elevs = atan2(zen, sqrt(north * north + w * w));
  *azs = atan2(-w, north);
  if (*azs < 0) *azs = *azs + TWOPI;
}

int moon(double time, double mchk) {
  double az, el, lat, lon, ra, dec;
  if (mchk == 0.0) return 1;
  lon = 116.5 * PI / 180.0;  // Boolardy
  lat = -26.7 * PI / 180.0;  // Boolardy
  moonradec(time, &ra, &dec, lat, lon);
  radec_azel(gst(time) - ra + lon, dec, lat, &az, &el);
  printf("moon az %f el %f\n", az * 180.0 / PI, el * 180.0 / PI);
  if (el < 0 && mchk == 1) return 1;
  if (el >= 0 && mchk == 2) return 1;
  return 0;
}

int sun(double time, double schk) {
  double az, el, lat, lon, ra, dec;
  if (schk == 0.0) return 1;
  lon = 116.5 * PI / 180.0;  // Boolardy
  lat = -26.7 * PI / 180.0;  // Boolardy
  sunradec(time, &ra, &dec);
  radec_azel(gst(time) - ra + lon, dec, lat, &az, &el);
  printf("sun az %f el %f\n", az * 180.0 / PI, el * 180.0 / PI);
  if (el < 0) return 1;
  return 0;
}

double fition(int np, int mode, double freqq[], double data[], double wt[], double dataout[], double *t150, double *spind, double *spcurv,
              double *ionabs, double *ionemm, double *speor, double *cov, double feor, double wid, int sig, int imode, double tt, double tau) {
  int k, i, j, nfit;
  double freq, spi, f150, f;
  static double wtt[32768];

  for (j = 0; j < np; j++)
    if (wt[j])
      wtt[j] = 1;
    else
      wtt[j] = 0;
  nfit = abs(mode);
  spi = 2.55;
  if ((freqq[0] + freqq[np - 1]) / 2.0 > 100.0)
    f150 = 150;
  else
    f150 = 75;
  for (i = 0; i < nfit; i++) bbrr[i] = 0;
  for (i = 0; i < nfit; i++) {
    for (j = 0; j < np; j++) {
      k = j + i * np;
      freq = freqq[j] / f150;
      if (i == 0) fitf[k] = pow(freq, -spi);
      if (i == 1) fitf[k] = pow(freq, -spi) * log(freq);
      if (i == 2 && nfit == 3) fitf[k] = pow(freq, -spi - 2.0) - 0 * 1.5 * pow(freq, -2.0);  // 1.5 = Telectron/Tsky
      if (i == 2 && nfit >= 4) fitf[k] = pow(freq, -spi - 2.0);
      if (i == 3) fitf[k] = pow(freq, -2.0);
      if (i == 4) fitf[k] = pow(freq, -spi) * (pow(log(freq), 2.0));
      //    if(i==4)fitf[k] = 1;   //  better for high band
      if (i == 5) fitf[k] = pow(freq, -spi) * (pow(log(freq), 3.0));
      
      // sig = 0 for B18
      if (sig >= 10) fitf[k] = pow(freq, -2.5 + i);                              // not quite as good for updn  best for dn only
      if (sig >= 20 && sig < 30) fitf[k] = pow(freq, -spi) * pow(log(freq), i);  // new quasi physical
      if (sig >= 30) fitf[k] = pow(log(freq), i);                                // new loglog
      // imode=0 for B18 RMS filter.
      if (imode) {
        if (i == 0) fitf[k] = pow(freq, -spi - 2.0);
        if (i == 1) fitf[k] = pow(freq, -2.0);
        if (i == 2) fitf[k] = 1;
      }
      //    fitf[k] = pow(freq,-spi + i);  // not quite as good for updn  best for dn only
      //    if((nfit>7 && mode>0) || (nfit>8 && mode<0))  fitf[k] = pow(freq,-spi-2 + i);  // not quite as good for updn  best for dn only

      if (mode < 0 && i && i == nfit - 1 && (sig % 10) == 0) {
        if (tau == 0)
          fitf[k] = -exp(-0.693 * (freqq[j] - feor) * (freqq[j] - feor) / (wid * wid * 0.25));  // negative for absorption
        else
          fitf[k] = -(1 - exp(-tau * exp(-(freqq[j] - feor) * (freqq[j] - feor) * (-log(-log((1 + exp(-tau)) / 2) / tau)) / (wid * wid * 0.25)))) /
                    (1 - exp(-tau));  // true width
      }
      if (mode < 0 && i && i == nfit - 1 && (sig % 10) == 1)
        fitf[k] = 0.5 * (tanh((1420 / freqq[j] - 1420 / feor) / (wid * 1420 / (feor * feor))) + 1.0);
      if (i == nfit - 1 && sig == 30) fitf[k] = fitf[k] / pow(freq, -spi);
      //    if(mode < 0 && i && i==nfit-1 && (sig%10) == 1) fitf[k] = 0.5*(tanh((1420/freqq[j]-1420/feor)/(wid*1420/(feor*feor-wid*wid/4.0))))+1.0);
      //    if(mode < 0 && i && i==nfit-1 && sig == 1) fitf[k] =
      //    sqrt((1420/freqq[j])/10.0)*0.5*(tanh((1420/freqq[j]-1420/feor)/(wid*1420/(feor*feor)))+1.0);
    }
  }
  if (sig >= 30)
    for (j = 0; j < np; j++) data[j] = log(data[j]);
  if (sig >= 30)
    for (j = 0; j < np; j++)
      if (wtt[j] > 0) wtt[j] = data[j] * data[j];  // to make equivalent to fitting without using non linear log(u+a)~log(u)+a/u
  polyfitr(nfit, np, data, mcalc, wtt, dataout, fitf, aarr, bbrr);
  if (sig >= 30)
    for (j = 0; j < np; j++)
      if (wtt[j] > 0) wtt[j] = 1;
  if (sig >= 30)
    for (j = 0; j < np; j++) {
      data[j] = exp(data[j]);
      dataout[j] = exp(dataout[j]);
    }
  *t150 = 0;
  for (i = 0; i < nfit; i++) {
    freq = 1;
    if (sig < 10) {
      if (i == 0) *t150 += bbrr[i] * pow(freq, -spi);
      if (i == 1) *t150 += bbrr[i] * pow(freq, -spi) * log(freq);
      if (i == 4) *t150 += bbrr[i] * pow(freq, -spi) * (pow(log(freq), 2.0));
      if (i == 2) *t150 += bbrr[i] * pow(freq, -spi - 2.0);
      if (i == 3) *t150 += bbrr[i] * pow(freq, -2.0);
      if (i == 5) *t150 += bbrr[i] * pow(freq, -spi) * (pow(log(freq), 3.0));
    }
    if (sig >= 10 && sig < 20) *t150 += bbrr[i] * pow(freq, -spi + i);
    if (sig >= 20 && sig < 30) *t150 += bbrr[i] * pow(freq, -spi) * pow(log(freq), i);
  }
  if (sig >= 30 && sig < 40) *t150 = exp(bbrr[0]);
  if (sig >= 30 && sig < 40 && mode < 1) {
    f = (exp(bbrr[0] + bbrr[1] * log(feor / f150))) / pow(feor / f150, -spi);
    bbrr[nfit - 1] = f * bbrr[nfit - 1];
    aarr[(nfit - 1) + (nfit - 1) * nfit] = sqrt(f) * aarr[(nfit - 1) + (nfit - 1) * nfit];
  }  // approx.
     //     *t150=bbrr[0];
  *spind = -spi + bbrr[1] / bbrr[0];
  *ionabs = -1e2 * bbrr[2] / bbrr[0];
  *ionemm = -bbrr[3] * bbrr[0] / bbrr[2];
  *spcurv = bbrr[4] / bbrr[0];
  *speor = bbrr[nfit - 1];
  *cov = aarr[(nfit - 1) + (nfit - 1) * nfit];
  if (nfit < 3 || (mode < 0 && nfit < 5)) *ionemm = 0;
  if (imode) {
    *spind = *spcurv = *speor = 0;
    *cov = 1;
    *ionabs = -1e2 * bbrr[0] / tt;
    *ionemm = -tt * bbrr[1] / bbrr[0];
  }  // tt comes from t150dn
     //  printf("cov %e nfit %d feor %f mode %d np %d wid %f sig %d wt %f\n",*cov,nfit,feor,mode,np,wid,sig,wt[np/2]);
     //  printf("b0 %e b1 %e mode %d b4 %f data %e %e %e\n",*b0,*b1,mode,*b4,data[np/2],wt[np/2],dataout[np/2]);
  return f150;
}

void eorsnr(int np, double freqq[], double data[], double wt[], double dataout[], double fstart, double fstop, double wid, int sig, int nfit,
            double seor, double tau, int mode, double fsmax, double fsmin, char info4[]) {
  double feor, t150, spind, spcurv, ionabs, ionemm, speor, cov, rms, snr, widd, widmax, widmaxx, snrmax, snrmaxx, maxfreq, sigmax, sigmaxx, rmsin,
      rmss, rmsss, wwid;
  double wtt[NDATA], snrr[NDATA], dataout1[NDATA], dataout2[NDATA], datasig[NDATA];
  int i, j;
  char info[256];
  fition(np, nfit, freqq, data, wt, dataout2, &t150, &spind, &spcurv, &ionabs, &ionemm, &speor, &cov, 150, wid, sig, 0, 0, tau);
  rmsin = rmscalc(np, data, dataout, wt);
  for (i = 0; i < np; i++) dataout1[i] = data[i] - dataout[i];
  if (fstop > 150)
    wwid = 60;
  else
    wwid = 30;
  snrmaxx = widmaxx = maxfreq = wid = snrmax = widmax = sigmax = sigmaxx = rmss = rmsss = 0;
  for (i = 0; i < np; i++) {
    feor = freqq[i];
    for (j = 0; j < np; j++) wtt[j] = wt[j];
    snrmax = widmax = sigmax = 0;
    for (widd = 10; widd <= wwid; widd += 0.1) {
      if (feor >= fstart + widd / 2.0 && feor <= fstop - widd / 2.0 && feor >= fsmin && feor <= fsmax) {
        fition(np, -(nfit + 1), freqq, data, wtt, dataout, &t150, &spind, &spcurv, &ionabs, &ionemm, &speor, &cov, feor, widd, sig, 0, 0, tau);
        rms = rmscalc(np, data, dataout, wtt);
        snr = speor / (sqrt(cov) * rms);
        if (seor == -2) {
          if (speor > 0)
            snr = 1.0 / rms;
          else
            snr = 0;
        }  // best fit with correct sign
        if (wtt[i] && snr > snrmax) {
          snrmax = snr;
          widmax = widd;
          sigmax = speor;
          rmss = rms;
        }
      }
    }
    if (snrmax > snrmaxx) {
      snrmaxx = snrmax;
      widmaxx = widmax;
      sigmaxx = sigmax;
      maxfreq = feor;
      rmsss = rmss;
    }
    snrr[i] = snrmax;
    printf("freq %f snrmax %f sigmax %f widmax %f fstart %f fstop %f\n", feor, snrmax, sigmax, widmax, fstart, fstop);
  }
  fition(np, -(nfit + 1), freqq, data, wtt, dataout, &t150, &spind, &spcurv, &ionabs, &ionemm, &speor, &cov, maxfreq, widmaxx, sig, 0, 0, tau);
  if (seor == -2) snrmaxx = speor / (sqrt(cov) * rmsss);
  printf("maxfreq %5.1f snrmax %5.1f sigmax %5.2f widmax %5.1f rmsin %6.4f rms %6.4f fstart %5.1f fstop %5.1f\n", maxfreq, snrmaxx, sigmaxx, widmaxx,
         rmsin, rmsss, fstart, fstop);
  sprintf(info, "freq %4.1f snr %4.1f sig %4.2f wid %4.2f tau %2.0f rmsin %6.4f rms %6.4f %2.0f - %2.0f\n", maxfreq, snrmaxx, sigmaxx, widmaxx, tau,
          rmsin, rmsss, fstart, fstop);
  // for(j=0;j<np;j++) dataout2[j] = data[j] - dataout2[j];
  for (j = 0; j < np; j++) dataout2[j] = data[j] - dataout[j];
  if (mode) {
    for (j = 0; j < np; j++) {
      if (tau)
        datasig[j] =
            -sigmaxx *
            (1 -
             exp(-tau * exp(-(freqq[j] - maxfreq) * (freqq[j] - maxfreq) * (-log(-log((1 + exp(-tau)) / 2) / tau)) / (widmaxx * widmaxx * 0.25)))) /
            (1 - exp(-tau));
      else {
        if ((sig % 10) == 0) datasig[j] = -sigmaxx * exp(-0.693 * (freqq[j] - maxfreq) * (freqq[j] - maxfreq) / (widmaxx * widmaxx * 0.25));
        if ((sig % 10) == 1) datasig[j] = sigmaxx * 0.5 * (tanh((1420 / freqq[j] - 1420 / maxfreq) / (widmaxx * 1420 / (maxfreq * maxfreq))) + 1.0);
      }
    }
    zplott2(np, freqq, dataout1, dataout2, datasig, snrr, wt, nfit, info, info4);
  } else
    zplott(np, freqq, dataout1, dataout2, snrr, wt, nfit, info);
}

double polyfitr(int npoly, int nfreq, double ddata[], double mcalc[], double wtt[], double dataout[], double fitfn[], long double aarr[],
                double bbrr[]) {
  int i, j, k, kk, m1, m2;
  double re, dd;
  for (i = 0; i < nfreq; i++) {
    kk = i * npoly;
    for (j = 0; j < npoly; j++) {
      mcalc[kk] = fitfun(j, i, nfreq, fitfn);
      kk++;
    }
  }
  for (i = 0; i < npoly; i++) {
    re = 0.0;
    m1 = i;
    for (k = 0; k < nfreq; k++) {
      dd = ddata[k];
      if (wtt[k] > 0) re += mcalc[m1] * dd * wtt[k];
      m1 += npoly;
    }
    bbrr[i] = re;
    for (j = 0; j < npoly; j++) {
      re = 0.0;
      m1 = i;
      m2 = j;
      for (k = 0; k < nfreq; k++) {
        if (wtt[k] > 0) re += mcalc[m1] * mcalc[m2] * wtt[k];
        m1 += npoly;
        m2 += npoly;
      }
      k = j + i * npoly;

      aarr[k] = re;
    }
  }
  qrd(aarr, npoly, bbrr);
  for (i = 0; i < nfreq; i++) {
    re = 0.0;
    for (j = 0; j < npoly; j++) { re += bbrr[j] * fitfun(j, i, nfreq, fitfn); }
    dd = re;
    dataout[i] = dd;
  }
  return bbrr[0];
}

double rmscalc(int n, double data[], double md[], double wt[]) {
  int i;
  double rms, m;
  rms = m = 0;
  for (i = 0; i < n; i++) {
    rms += wt[i] * (data[i] - md[i]) * (data[i] - md[i]);
    m += wt[i];
  }
  return sqrt(rms / m);
}

int chisolve(int np, double freqq[], double data[], double wt[], double dataout[], int sig, int nfit, double fstart, double fstop, double fsmax,
             double fsmin) {
  double tau, aasig, ffsig, wwidsig, ttau, rmsfit, dfreq, damp, dwid, dtau, snr;
  int phy, n, i, j, k, l, m, min;
  static unsigned char map[65536];  // 16x16x16x16
  char rslts[256];
  for (i = 0; i < 65536; i++) map[i] = 255;
  tau = 7;
  dfreq = 0.1;
  damp = 0.02;
  dwid = 0.2;
  dtau = 0.5;
  if (sig == 0)
    phy = 2;
  else
    phy = 0;
  snr = chifit(0, np, freqq, data, wt, dataout, nfit, tau, &aasig, &ffsig, &wwidsig, &ttau, &rmsfit, phy, fsmax, fsmin, map, dfreq, damp, dwid, dtau);
  printf("aasig %f ffsig %f wwidsig %f ttau %f\n", aasig, ffsig, wwidsig, ttau);
  for (n = 1; n <= 4; n++)
    chifit(n, np, freqq, data, wt, dataout, nfit, tau, &aasig, &ffsig, &wwidsig, &ttau, &rmsfit, phy, fsmax, fsmin, map, dfreq, damp, dwid, dtau);

  printf("amp down\n");
  for (j = 0; j < 16; j++) {
    for (i = 0; i < 16; i++) {
      min = 99;
      {
        for (k = 0; k < 16; k++)
          for (l = 0; l < 16; l++) {
            m = i + j * 16 + k * 16 * 16 + l * 16 * 16 * 16;
            if (map[m] < min) min = map[m];
          }
      }
      printf("%d", min / 10);
    }
    printf("\n");
  }
  printf("frequency across\n");

  printf("1\n");

  printf("width down\n");
  for (k = 0; k < 16; k++) {
    for (i = 0; i < 16; i++) {
      min = 99;
      {
        for (j = 0; j < 16; j++)
          for (l = 0; l < 16; l++) {
            m = i + j * 16 + k * 16 * 16 + l * 16 * 16 * 16;
            if (map[m] < min) min = map[m];
          }
      }
      printf("%d", min / 10);
    }
    printf("\n");
  }
  printf("frequency across\n");

  printf("2\n");

  printf("tau down\n");
  for (i = 0; l < 16; l++) {
    for (i = 0; i < 16; i++) {
      min = 99;
      {
        for (j = 0; j < 16; j++)
          for (k = 0; k < 16; k++) {
            m = i + j * 16 + k * 16 * 16 + l * 16 * 16 * 16;
            if (map[m] < min) min = map[m];
          }
      }
      printf("%d", min / 10);
    }
    printf("\n");
  }
  printf("frequency across\n");

  printf("3\n");

  printf("width down\n");
  for (k = 0; k < 16; k++) {
    for (j = 0; j < 16; j++) {
      min = 99;
      {
        for (i = 0; i < 16; i++)
          for (l = 0; l < 16; l++) {
            m = i + j * 16 + k * 16 * 16 + l * 16 * 16 * 16;
            if (map[m] < min) min = map[m];
          }
      }
      printf("%d", min / 10);
    }
    printf("\n");
  }
  printf("amplitude across\n");

  printf("4\n");

  printf("tau down\n");
  for (l = 0; l < 16; l++) {
    for (j = 0; j < 16; j++) {
      min = 99;
      {
        for (i = 0; i < 16; i++)
          for (k = 0; k < 16; k++) {
            m = i + j * 16 + k * 16 * 16 + l * 16 * 16 * 16;
            if (map[m] < min) min = map[m];
          }
      }
      printf("%d", min / 10);
    }
    printf("\n");
  }
  printf("amplitude across\n");

  printf("5\n");

  printf("tau down\n");
  for (l = 0; l < 16; l++) {
    for (k = 0; k < 16; k++) {
      min = 99;
      {
        for (i = 0; i < 16; i++)
          for (j = 0; j < 16; j++) {
            m = i + j * 16 + k * 16 * 16 + l * 16 * 16 * 16;
            if (map[m] < min) min = map[m];
          }
      }
      printf("%d", min / 10);
    }
    printf("\n");
  }
  printf("width across\n");

  printf("6\n");
  sprintf(rslts, "fcen %5.2f amp %5.2f wid %5.2f tau %5.2f nterms %d phy %d fstart %5.2f fstop %5.2f fsmin %5.2f fsmax %5.2f snr %5.2f", ffsig, aasig,
          wwidsig, ttau, nfit, phy, fstart, fstop, fsmin, fsmax, snr);
  plotq(map, ffsig, dfreq, aasig, damp, wwidsig, dwid, ttau, dtau, rslts);
  return 0;
}

double chifit(int mode, int np, double freqq[], double data[], double wt[], double dataout[], int nfit, double tau, double *aasig, double *ffsig,
              double *wwidsig, double *ttau, double *rmsfit, int phy, double fsmax, double fsmin, unsigned char map[], double dfreq, double damp,
              double dwid, double dtau) {
  double feor, speor, cov, rms, snr, widd, widmax, widmaxx, snrmax, snrmaxx, maxfreq, sigmax, sigmaxx, rmsin, rmss, rmsss, tttau, tau1, tau2, taumax,
      taumaxx;
  double wtt[NDATA], data1[NDATA], sig, chi, chi0, nn;
  int i, j, k;
  rmss = 0;
  tau1 = 1;
  tau2 = tau + 5;
  nn = 0;
  for (i = 0; i < np; i++) nn += wt[i];
  if (mode == 0) {
    fit(np, phy, nfit, freqq, data, wt, dataout, &speor, &cov, 0, 0, 0);
    rmsin = rmscalc(np, data, dataout, wt);
    snrmaxx = widmaxx = maxfreq = snrmax = widmax = sigmax = sigmaxx = rmss = rmsss = taumax = taumaxx = 0;
    for (feor = fsmin; feor <= fsmax; feor += dfreq) {
      for (j = 0; j < np; j++) wtt[j] = wt[j];
      snrmax = widmax = sigmax = 0;
      for (tttau = tau1; tttau <= tau2; tttau += dtau) {
        for (widd = 10; widd <= 20; widd += dwid) {
          fit(np, phy + 1, nfit, freqq, data, wtt, dataout, &speor, &cov, feor, widd, tttau);
          rms = rmscalc(np, data, dataout, wtt);
          //   snr = speor/(sqrt(cov)*rms);
          snr = 1.0 / rms;
          if (snr > snrmax) {
            snrmax = snr;
            widmax = widd;
            sigmax = speor;
            taumax = tttau;
            rmss = rms;
          }
        }
      }
      printf("freq %f wtt %f snrmax %f sigmax %f widmax %f taumax %f\n", feor, wtt[i], snrmax, sigmax, widmax, taumax);
      if (snrmax > snrmaxx) {
        snrmaxx = snrmax;
        widmaxx = widmax;
        sigmaxx = sigmax;
        maxfreq = feor;
        taumaxx = taumax;
        rmsss = rmss;
      }
    }
    fit(np, phy + 1, nfit, freqq, data, wtt, dataout, &speor, &cov, maxfreq, widmaxx, taumaxx);
    rms = rmscalc(np, data, dataout, wtt);
    snr = speor / (sqrt(cov) * rms);
    printf("maxfreq %5.1f snrmax %5.1f sigmax %5.2f widmax %5.1f taumax %5.1f rmsin %6.4f rms %6.4f snr %5.1f\n", maxfreq, snrmaxx, sigmaxx, widmaxx,
           taumaxx, rmsin, rmsss, snr);
    *aasig = sigmaxx;
    *ffsig = maxfreq;
    *wwidsig = widmaxx;
    *ttau = taumaxx;
    *rmsfit = rmsss;

    return snr;
  }
  if (mode == 1) {
    printf("mode 1\n");
    for (sig = *aasig - 0.2; sig <= *aasig + 0.2; sig += damp) {  // amp error
      snrmax = 0;
      for (feor = *ffsig - 1; feor <= *ffsig + 1; feor += dfreq) {
        for (j = 0; j < np; j++) wtt[j] = wt[j];
        for (tttau = tau1; tttau <= tau2; tttau += dtau) {
          for (widd = *wwidsig - 2; widd <= *wwidsig + 2; widd += 0.2) {
            for (k = 0; k < np; k++)
              data1[k] = data[k] + sig *
                                       (1 - exp(-tttau * exp(-(freqq[k] - feor) * (freqq[k] - feor) * (-log(-log((1 + exp(-tttau)) / 2) / tttau)) /
                                                             (widd * widd * 0.25)))) /
                                       (1 - exp(-tttau));
            fit(np, phy, nfit, freqq, data1, wtt, dataout, &speor, &cov, feor, widd, tttau);
            rms = rmscalc(np, data1, dataout, wtt);
            mapp(map, feor, *ffsig, dfreq, sig, *aasig, damp, widd, *wwidsig, dwid, tttau, *ttau, dtau,
                 nn * ((rms / (*rmsfit)) * (rms / (*rmsfit)) - 1));
            snr = 1.0 / rms;
            if (snr > snrmax) {
              snrmax = snr;
              rmss = rms;
            }
          }
        }
      }
      chi0 = 1;
      chi = rmss / (*rmsfit);
      chi = nn * (chi * chi - 1);
      printf("sig %f rmss %f chi0 %f chi %f\n", sig, rmss, chi0, chi);
    }
  }

  if (mode == 2) {  // wid error
    printf("mode 2\n");
    for (widd = *wwidsig - 2; widd <= *wwidsig + 2; widd += dwid) {
      snrmax = 0;
      for (feor = *ffsig - 1; feor <= *ffsig + 1; feor += dfreq) {
        for (j = 0; j < np; j++) wtt[j] = wt[j];
        for (tttau = tau1; tttau <= tau2; tttau += dtau) {
          fit(np, phy + 1, nfit, freqq, data, wtt, dataout, &speor, &cov, feor, widd, tttau);
          rms = rmscalc(np, data, dataout, wtt);
          mapp(map, feor, *ffsig, dfreq, speor, *aasig, damp, widd, *wwidsig, dwid, tttau, *ttau, dtau,
               nn * ((rms / (*rmsfit)) * (rms / (*rmsfit)) - 1));
          snr = 1.0 / rms;
          if (snr > snrmax) {
            snrmax = snr;
            rmss = rms;
          }
        }
      }
      chi0 = 1;
      chi = rmss / (*rmsfit);
      chi = nn * (chi * chi - 1);
      printf("widd %f rmss %f chi0 %f chi %f\n", widd, rmss, chi0, chi);
    }
  }

  if (mode == 3) {  // feor error
    printf("mode 3\n");
    for (feor = *ffsig - 1; feor <= *ffsig + 1; feor += dfreq) {
      snrmax = 0;
      for (j = 0; j < np; j++) wtt[j] = wt[j];
      for (tttau = tau1; tttau <= tau2; tttau += dtau) {
        for (widd = *wwidsig - 2; widd <= *wwidsig + 2; widd += dwid) {
          fit(np, phy + 1, nfit, freqq, data, wtt, dataout, &speor, &cov, feor, widd, tttau);
          rms = rmscalc(np, data, dataout, wtt);
          mapp(map, feor, *ffsig, dfreq, speor, *aasig, damp, widd, *wwidsig, dwid, tttau, *ttau, dtau,
               nn * ((rms / (*rmsfit)) * (rms / (*rmsfit)) - 1));
          snr = 1.0 / rms;
          if (snr > snrmax) {
            snrmax = snr;
            rmss = rms;
          }
        }
      }
      chi0 = 1;
      chi = rmss / (*rmsfit);
      chi = nn * (chi * chi - 1);
      printf("feor %f rmss %f chi0 %f chi %f\n", feor, rmss, chi0, chi);
    }
  }

  if (mode == 4) {  // tau error
    printf("mode 4\n");
    for (tttau = tau1; tttau <= tau2; tttau += dtau) {
      snrmax = 0;
      for (j = 0; j < np; j++) wtt[j] = wt[j];
      for (feor = *ffsig - 1; feor <= *ffsig + 1; feor += dfreq) {
        for (widd = *wwidsig - 2; widd <= *wwidsig + 2; widd += dwid) {
          fit(np, phy + 1, nfit, freqq, data, wtt, dataout, &speor, &cov, feor, widd, tttau);
          rms = rmscalc(np, data, dataout, wtt);
          mapp(map, feor, *ffsig, dfreq, speor, *aasig, damp, widd, *wwidsig, dwid, tttau, *ttau, dtau,
               nn * ((rms / (*rmsfit)) * (rms / (*rmsfit)) - 1));
          snr = 1.0 / rms;
          if (snr > snrmax) {
            snrmax = snr;
            rmss = rms;
          }
        }
      }
      chi0 = 1;
      chi = rmss / (*rmsfit);
      chi = nn * (chi * chi - 1);
      printf("tau %f rmss %f chi0 %f chi %f\n", tttau, rmss, chi0, chi);
    }
  }
  return 0;
}

void mapp(unsigned char map[], double feor, double feor0, double dfeor, double sig, double sig0, double dsig, double wid, double wid0, double dwid,
          double tau, double tau0, double dtau, double chi) {
  int i, j, k, l, m, n;
  i = (feor - feor0) / dfeor + 8.5;
  if (i < 0) i = 0;
  if (i > 15) i = 15;
  j = (sig - sig0) / dsig + 8.5;
  if (j < 0) j = 0;
  if (j > 15) j = 15;
  k = (wid - wid0) / dwid + 8.5;
  if (k < 0) k = 0;
  if (k > 15) k = 15;
  l = (tau - tau0) / dtau + 8.5;
  if (l < 0) l = 0;
  if (l > 15) l = 15;
  m = i + j * 16 + k * 16 * 16 + l * 16 * 16 * 16;
  n = 10.0 * chi;
  if (n < 0) n = 0;
  if (n > 255) n = 255;
  map[m] = n;
}

void fit(int np, int mode, int nfit, double freqq[], double data[], double wt[], double dataout[], double *speor, double *cov, double feor,
         double wid, double tau) {
  int k, i, j, mfit;
  double freq;
  static double wtt[32768];

  for (j = 0; j < np; j++)
    if (wt[j])
      wtt[j] = 1;
    else
      wtt[j] = 0;
  mfit = nfit;
  if (mode & 1) mfit = nfit + 1;
  for (i = 0; i < 8; i++) bbrr[i] = 0;
  for (i = 0; i < mfit; i++) {
    for (j = 0; j < np; j++) {
      k = j + i * np;
      freq = freqq[j] / 75;
      fitf[k] = pow(freq, -2.5 + i);
      if (mode & 2) {
        if (i == 0) fitf[k] = pow(freq, -2.55);
        if (i == 1) fitf[k] = pow(freq, -2.55) * log(freq);
        if (i == 2 && nfit == 3) fitf[k] = pow(freq, -2.55 - 2.0);
        if (i == 2 && nfit >= 4) fitf[k] = pow(freq, -2.55 - 2.0);
        if (i == 3) fitf[k] = pow(freq, -2.0);
        if (i == 4) fitf[k] = pow(freq, -2.55) * (pow(log(freq), 2.0));
      }

      if ((mode & 1) && i == mfit - 1)
        fitf[k] = -(1 - exp(-tau * exp(-(freqq[j] - feor) * (freqq[j] - feor) * (-log(-log((1 + exp(-tau)) / 2) / tau)) / (wid * wid * 0.25)))) /
                  (1 - exp(-tau));  // true width
    }
  }
  polyfitr(mfit, np, data, mcalc, wtt, dataout, fitf, aarr, bbrr);
  *speor = bbrr[mfit - 1];
  *cov = aarr[(mfit - 1) + (mfit - 1) * mfit];
}

double fitfun(int j, int i, int nfreq, double fitfn[]) { return fitfn[i + j * nfreq]; }

void qrd(long double a[], int n, double b[]) {
  int i, j, k;
  int pi, pk;
  long double c[100], d[100];
  long double scale, sigma, sum, tau;
  static long double qt[30][30], u[30][30];  // aa[30][30];

  //  for(i=0;i<n;i++)for(j=0;j<n;j++) aa[i][j]=a[j+i*n]; // save for check

  for (k = 0; k < n - 1; k++) {
    scale = 0.0;
    for (i = k; i < n; i++)
      if (fabsl(a[k + i * n])) scale = fabsl(a[k + i * n]);
    if (scale == 0.0) {
      printf("singular1\n");
      c[k] = d[k] = 0.0;
    } else {
      for (i = k; i < n; i++) a[k + i * n] /= scale;
      for (sum = 0.0, i = k; i < n; i++) sum += a[k + i * n] * a[k + i * n];
      if (a[k + k * n] > 0)
        sigma = sqrtl(sum);
      else
        sigma = -sqrtl(sum);
      a[k + k * n] += sigma;
      c[k] = sigma * a[k + k * n];
      d[k] = -scale * sigma;
      for (j = k + 1; j < n; j++) {
        for (sum = 0.0, i = k; i < n; i++) sum += a[k + i * n] * a[j + i * n];
        tau = sum / c[k];
        for (i = k; i < n; i++) a[j + i * n] -= tau * a[k + i * n];
      }
    }
  }
  d[n - 1] = a[n - 1 + (n - 1) * n];
  if (d[n - 1] == 0.0) printf("singular2\n");
  for (i = 0; i < n; i++) {  // Form QT explicitly.
    for (j = 0; j < n; j++) qt[i][j] = 0.0;
    qt[i][i] = 1.0;
  }
  for (k = 0; k < n - 1; k++) {
    if (c[k] != 0.0) {
      for (j = 0; j < n; j++) {
        sum = 0.0;
        for (i = k; i < n; i++) sum += a[k + i * n] * qt[i][j];
        sum /= c[k];
        for (i = k; i < n; i++) qt[i][j] -= sum * a[k + i * n];
      }
    }
  }
  for (j = 0; j < n - 1; j++) {
    for (sum = 0, i = j; i < n; i++) sum += a[j + i * n] * b[i];
    tau = sum / c[j];
    for (i = j; i < n; i++) b[i] -= tau * a[j + i * n];
  }

  b[n - 1] /= d[n - 1];
  for (i = n - 2; i >= 0; i--) {
    for (sum = 0, j = i + 1; j < n; j++) sum += a[j + i * n] * b[j];
    b[i] = (b[i] - sum) / d[i];
  }

  // int Upper_Triangular_Inverse(double *U, int n)

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (j > i)
        u[i][j] = a[j + i * n];
      else
        u[i][j] = 0;
      if (i == j) u[i][i] = d[i];
    }
  }

  for (k = 0; k < n; k++) {
    if (u[k][k] == 0)
      return;
    else
      u[k][k] = 1.0 / u[k][k];
  }

  for (i = n - 2, pi = (n - 2); i >= 0; pi -= 1, i--) {
    for (j = n - 1; j > i; j--) {
      sum = 0.0;
      for (k = i + 1, pk = pi + 1; k <= j; pk += 1, k++) { sum += u[pi][k] * u[pk][j]; }
      u[pi][j] = -u[pi][i] * sum;
    }
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      sum = 0;
      for (k = 0; k < n; k++) sum += u[i][k] * qt[k][j];
      a[j + i * n] = sum;
    }
  }
  /*
   for(i=0;i<n;i++){
    printf("sum ");
    for(j=0;j<n;j++) {
     sum=0;
    for(k=0;k<n;k++) sum +=a[k+i*n]*aa[k][j];
        printf("%10.3Le ",sum);
     }
     printf("\n");
    }
  */
  // A^{-1} = R^{-1}Q^^{-1} = R^{-1}QT   - see Press 3rd edition
}

void outspec(int nn[], double freqq[], double data[], double data2[], double data3[], double wt[], int nplot, int out, char title[])
//   writeout spectrum
{
  char txt[256];
  int k, np;
  FILE *file;

  sprintf(txt, "spe0.txt");
  if ((file = fopen(txt, "w")) == NULL) {
    printf("cannot open %s:\n", txt);
    return;
  }
  np = nn[nplot - 1];
  for (k = 0; k < np; k++) {
    if (out == 1)
      fprintf(file, "freq %8.4f resid %8.4f K wt %8.4f\n", freqq[k + (nplot - 1) * NDATA], data[k + (nplot - 1) * NDATA],
              wt[k + (nplot - 1) * NDATA]);
    if (out == 2)
      fprintf(file, "freq %10.6f tantenna %10.6f K skymodel %10.6f K wt %1.0f %s\n", freqq[k + (nplot - 1) * NDATA], data[k + (nplot - 1) * NDATA],
              0.0, wt[k + (nplot - 1) * NDATA], title);
    if (out == 3)
      fprintf(file, "freq %10.6f tantenna %10.6f K skymodel %10.6f K resid %10.6f K wt %1.0f %s\n", freqq[k + (nplot - 1) * NDATA], data2[k],
              data3[k], data[k + (nplot - 1) * NDATA], wt[k + (nplot - 1) * NDATA], title);
    if (out == 4)
      fprintf(file, "freq %10.6f tantenna %10.6f K skymodel %10.6f K wt %1.0f %s\n", freqq[k + (nplot - 1) * NDATA], data2[k], 0.0,
              wt[k + (nplot - 1) * NDATA], title);
  }
  fclose(file);
}

double gauss(void) {
  double v1, v2, r, fac, aamp, vv1;
  v1 = r = 0.0;
  while (r > 1.0 || r == 0.0) {
    v1 = 2.0 * (rand() / 2147483648.0) - 1.0;
    v2 = 2.0 * (rand() / 2147483648.0) - 1.0;
    r = v1 * v1 + v2 * v2;
  }
  fac = sqrt(-2.0 * log(r) / r);
  vv1 = v1 * fac;
  aamp = vv1;
  return (aamp);
}

void zplott(int np, double freqq[], double datain[], double data[], double snr[], double wt[], int nfit, char info[]) {
  char txt[256];
  int k, m, iter;
  double x, y, xp, yp, dmax, dmin, f, totpp, scale, step, h, s, b, z;
  double xoffset, yoffset, fstart, fstop, rmsin, rms;
  FILE *file, *file1;

  sprintf(txt, "spez.pos");
  if ((file = fopen(txt, "w")) == NULL) {
    printf("cannot open %s:\n", txt);
    return;
  }

  sprintf(txt, "spez.txt");
  if ((file1 = fopen(txt, "w")) == NULL) {
    printf("cannot open %s:\n", txt);
    return;
  }
  for (k = 0; k < np; k++) {
    if (wt[k])
      m = 1;
    else
      m = 0;
    fprintf(file1, "freq %8.5f MHz spect %8.1f mK wt %d snr %8.5f MHz\n", freqq[k], data[k] * 1e3, m, snr[k]);
  }
  fclose(file1);

  fprintf(file, "%%!PS-Adobe-3.0 EPSF-3.0\n%c%cBoundingBox:  0 0 612 700\n%c%cEndProlog\n", '%', '%', '%', '%');
  fprintf(file, "1 setlinewidth\n");
  fprintf(file, "/Times-Roman findfont\n 12 scalefont\n setfont\n");

  rms = m = 0;
  for (k = 0; k < np; k++) {
    rms += wt[k] * data[k] * data[k];
    m += wt[k];
  }
  rms = sqrt(rms / m);
  rmsin = m = 0;
  for (k = 0; k < np; k++) {
    rmsin += wt[k] * datain[k] * datain[k];
    m += wt[k];
  }
  rmsin = sqrt(rmsin / m);

  fstart = fstop = 0;
  xoffset = 80.0;
  xoffset = 60.0;
  yoffset = 100.0;
  if (freqq[np / 2] > 100) {
    dmin = -100e-3;
    dmax = 100e-3;
  } else {
    dmax = 200e-3;
    dmin = -200e-3;
  }
  scale = dmax - dmin;
  for (y = 0; y < 2; y++)
    fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset, y * 480 + yoffset, xoffset + 400.0,
            y * 480 + yoffset);
  fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset, yoffset, xoffset, 480.0 + yoffset);
  fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset + 400.0, yoffset, xoffset + 400.0,
          480.0 + yoffset);
  fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset + 0.0, yoffset, xoffset + 0.0,
          480.0 + yoffset);
  fstart = freqq[0];
  if ((int)fstart == 90) fstart = 100;
  for (iter = 0; iter < 3; iter++) {
    fstop = (int)(freqq[np - 1] + 0.5);  // set to nearest MHz
    yp = 0;
    xp = 0;
    totpp = 0;
    for (k = 0; k < np; k++) {
      h = s = b = 0;
      x = (freqq[k] - fstart) * 400.0 / (fstop - fstart);
      if (iter == 0) totpp = (data[k] - dmin) / (2 * scale);
      if (iter == 1) totpp = (datain[k] - dmin) / (2 * scale) + 0.5;
      y = totpp * 240;
      if (iter == 2) y = 240.0 + snr[k] * 240.0 / 50.0;
      if (y > 480.0) y = 480.0;
      if (y < 0.0) y = 0.0;
      if (k == 0 || x < xp) {
        xp = x;
        yp = y;
      }
      h = s = b = 0;
      if ((k && wt[k] > 0.0 && wt[k - 1] > 0.0) || (k && iter == 2)) {
        fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n %5.2f %5.2f %5.2f sethsbcolor stroke\n", xp + xoffset, yp + yoffset,
                x + xoffset, y + yoffset, h, s, b);
      }
      xp = x;
      yp = y;
    }
  }
  fprintf(file, "/Times-Roman findfont\n 14 scalefont\n setfont\n");
  step = 5.0;
  if (fstop - fstart > 60) step = 10;
  if (fstop - fstart > 100) step = 20;
  for (f = fstart; f <= fstop + step * 0.1; f += step) {
    x = (f - fstart) * 400.0 / (fstop - fstart);
    y = 0;
    fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, y + yoffset, x + xoffset,
            y + 5.0 + yoffset);
    //        sprintf(txt, "%5.1f", f);
    sprintf(txt, "%4.0f", f);
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 15.0, y - 20.0 + yoffset, txt);
  }

  for (z = 7; z <= 50; z += 1) {
    f = 1420.0 / (z + 1);
    x = (f - fstart) * 400.0 / (fstop - fstart);
    if (x >= 0 && x <= 400.0) {
      y = 0;
      fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, y + yoffset + 480, x + xoffset,
              y - 5.0 + yoffset + 480);
      sprintf(txt, "%4.0f", z);
      fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 15.0, y - 15.0 + yoffset + 480 + 20, txt);
    }
  }
  sprintf(txt, " Z ");
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 220.0, 70.0 + 530.0, txt);

  for (f = 0; f <= 50; f += 5) {
    x = 0;
    y = f / 50.0;
    y = y * 240.0 + 240.0;
    fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, y + yoffset, x + xoffset + 5,
            y + yoffset);
    sprintf(txt, "%4.0f", f);
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 25.0, y - 1.0 + yoffset, txt);
  }

  sprintf(txt, "SNR");
  fprintf(file, "%6.2f %6.2f moveto\n 90 rotate\n (%s) show\n -90 rotate\n", xoffset - 25.0, 150.0 + 300.0, txt);

  for (iter = 0; iter < 2; iter++) {
    for (f = dmin; f < dmax - 10e-3; f += dmax / 5) {
      x = 0;
      y = (f - dmin) * 0.5 / (dmax - dmin);
      y = y * 240.0 + iter * 120;
      fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, y + yoffset, x + xoffset + 5,
              y + yoffset);
      //        fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n",
      //           x + xoffset,y-240.0/2+yoffset,x + xoffset + 400,y-240.0/2+yoffset);

      {
        sprintf(txt, "%4.0f", f * 1e3);
        fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 28.0, y - 1.0 + yoffset, txt);
      }
    }
  }
  fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, 240.0 + yoffset, x + xoffset + 400,
          240.0 + yoffset);

  //    sprintf(txt, " Spectrum for 2015_204 to 2015_266");
  //    sprintf(txt, " Spectrum for ");
  //    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n",xoffset + 50.0 , yoffset + 35.0 + 120, txt);

  sprintf(txt, " Residual with %d terms removed rms = %3.0f mK", nfit, rmsin * 1e3);
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 30.0, yoffset + 20.0 + 70, txt);

  //    sprintf(txt, " Spectrum for ");
  //    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n",xoffset + 50.0 , yoffset + 35.0, txt);

  sprintf(txt, " Residual with %d terms plus signature removed rms = %3.0f mK", nfit, rms * 1e3);
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 30.0, yoffset + 20.0, txt);

  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 10.0, 45.0, info);

  sprintf(txt, " Temperature (mK)");
  //    fprintf(file, "%6.2f %6.2f moveto\n 90 rotate\n (%s) show\n -90 rotate\n",xoffset -52.0, 265.0, txt);
  fprintf(file, "%6.2f %6.2f moveto\n 90 rotate\n (%s) show\n -90 rotate\n", xoffset - 30.0, 150.0, txt);

  sprintf(txt, "Frequency (MHz)");
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 160.0, 65.0, txt);
  fprintf(file, "showpage\n%c%cTrailer\n", '%', '%');
  fclose(file);
}

void zplott2(int np, double freqq[], double datain[], double data[], double datasig[], double snr[], double wt[], int nfit, char info[],
             char info4[]) {
  char txt[256];
  int k, m, iter;
  double x, y, xp, yp, dmax, dmin, f, totpp, scale, step, h, s, b, z;
  double xoffset, yoffset, fstart, fstop, rmsin, rms;
  FILE *file, *file1;

  sprintf(txt, "spez.pos");
  if ((file = fopen(txt, "w")) == NULL) {
    printf("cannot open %s:\n", txt);
    return;
  }

  sprintf(txt, "spez.txt");
  if ((file1 = fopen(txt, "w")) == NULL) {
    printf("cannot open %s:\n", txt);
    return;
  }
  for (k = 0; k < np; k++) {
    if (wt[k])
      m = 1;
    else
      m = 0;
    fprintf(file1, "freq %8.5f MHz spect %8.1f mK wt %d snr %8.5f MHz\n", freqq[k], data[k] * 1e3, m, snr[k]);
  }
  fclose(file1);

  fprintf(file, "%%!PS-Adobe-3.0 EPSF-3.0\n%c%cBoundingBox:  0 0 612 700\n%c%cEndProlog\n", '%', '%', '%', '%');
  fprintf(file, "1 setlinewidth\n");
  fprintf(file, "/Times-Roman findfont\n 12 scalefont\n setfont\n");

  rms = m = 0;
  for (k = 0; k < np; k++) {
    rms += wt[k] * data[k] * data[k];
    m += wt[k];
  }
  rms = sqrt(rms / m);
  rmsin = m = 0;
  for (k = 0; k < np; k++) {
    rmsin += wt[k] * datain[k] * datain[k];
    m += wt[k];
  }
  rmsin = sqrt(rmsin / m);

  fstart = fstop = 0;
  xoffset = 80.0;
  xoffset = 60.0;
  yoffset = 100.0;
  if (freqq[np / 2] > 100 && freqq[0] > 80) {
    dmin = -100e-3;
    dmax = 100e-3;
  } else {
    dmax = 400e-3;
    dmin = -400e-3;
  }
  scale = dmax - dmin;
  for (y = 0; y < 2; y++)
    fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset, y * 480 + yoffset, xoffset + 400.0,
            y * 480 + yoffset);
  fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset, yoffset, xoffset, 480.0 + yoffset);
  fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset + 400.0, yoffset, xoffset + 400.0,
          480.0 + yoffset);
  fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset + 0.0, yoffset, xoffset + 0.0,
          480.0 + yoffset);
  fstart = freqq[0];
  //    if((int)fstart==90) fstart = 100;
  for (iter = 0; iter < 4; iter++) {
    fstop = (int)(freqq[np - 1] + 0.5);  // set to nearest MHz
    yp = 0;
    xp = 0;
    totpp = 0;
    for (k = 0; k < np; k++) {
      h = s = b = 0;
      x = (freqq[k] - fstart) * 400.0 / (fstop - fstart);
      if (iter == 0) totpp = (data[k] - dmin) / scale;
      if (iter == 1) totpp = (datain[k] - dmin) / scale;
      if (iter == 2) totpp = (datasig[k] - dmin) / scale;
      y = totpp * 80 + 80 * iter;
      if (iter == 3) y = 240.0 + snr[k] * 240.0 / 50.0;
      if (y > 480.0) y = 480.0;
      if (y < 0.0) y = 0.0;
      if (k == 0 || x < xp) {
        xp = x;
        yp = y;
      }
      h = s = b = 0;
      if ((k && wt[k] > 0.0 && wt[k - 1] > 0.0) || (k && iter == 3)) {
        fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n %5.2f %5.2f %5.2f sethsbcolor stroke\n", xp + xoffset, yp + yoffset,
                x + xoffset, y + yoffset, h, s, b);
      }
      xp = x;
      yp = y;
    }
  }
  fprintf(file, "/Times-Roman findfont\n 14 scalefont\n setfont\n");
  step = 5.0;
  if (fstop - fstart > 60) step = 10;
  if (fstop - fstart > 100) step = 20;
  for (f = fstart; f <= fstop + step * 0.1; f += step) {
    x = (f - fstart) * 400.0 / (fstop - fstart);
    y = 0;
    fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, y + yoffset, x + xoffset,
            y + 5.0 + yoffset);
    //        sprintf(txt, "%5.1f", f);
    sprintf(txt, "%4.0f", f);
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 15.0, y - 20.0 + yoffset, txt);
  }

  for (z = 7; z <= 50; z += 1) {
    f = 1420.0 / (z + 1);
    x = (f - fstart) * 400.0 / (fstop - fstart);
    if (x >= 0 && x <= 400.0) {
      y = 0;
      fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, y + yoffset + 480, x + xoffset,
              y - 5.0 + yoffset + 480);
      if (z <= 22 || ((int)z) % 2 == 0) {
        sprintf(txt, "%4.0f", z);
        fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 15.0, y - 15.0 + yoffset + 480 + 20, txt);
      }
    }
  }
  sprintf(txt, " Z ");
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 220.0, 70.0 + 530.0, txt);

  for (f = 0; f <= 50; f += 5) {
    x = 0;
    y = f / 50.0;
    y = y * 240.0 + 240.0;
    fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, y + yoffset, x + xoffset + 5,
            y + yoffset);
    sprintf(txt, "%4.0f", f);
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 25.0, y - 1.0 + yoffset, txt);
  }

  sprintf(txt, "SNR");
  fprintf(file, "%6.2f %6.2f moveto\n 90 rotate\n (%s) show\n -90 rotate\n", xoffset - 25.0, 150.0 + 300.0, txt);

  for (iter = 0; iter < 3; iter++) {
    for (f = dmin; f < dmax - 10e-3; f += dmax / 2) {
      x = 0;
      y = (f - dmin) / (dmax - dmin);
      y = y * 80.0 + iter * 80;
      fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, y + yoffset, x + xoffset + 5,
              y + yoffset);
      //        fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n",
      //           x + xoffset,y-240.0/2+yoffset,x + xoffset + 400,y-240.0/2+yoffset);

      if (iter < 2) {
        sprintf(txt, "%4.0f", f * 1e3);
        fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 28.0, y - 1.0 + yoffset, txt);
      } else {
        sprintf(txt, "%4.0f", f * 1e3 - 0 * dmax * 1e3);
        fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 28.0, y - 1.0 + yoffset, txt);
      }
    }
  }
  fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, 240.0 + yoffset, x + xoffset + 400,
          240.0 + yoffset);

  sprintf(txt, " Best fit signature");
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 30.0, yoffset + 0.0 + 160, txt);

  sprintf(txt, " Residual with %d terms removed rms = %3.0f mK", nfit, rmsin * 1e3);
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 30.0, yoffset + 0.0 + 80, txt);

  sprintf(txt, " Residual with %d terms plus signature removed rms = %3.0f mK", nfit, rms * 1e3);
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 30.0, yoffset + 20.0, txt);

  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 10.0, 45.0, info);

  if (strlen(info4) > 4) fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 250.0, 550.0, info4);

  sprintf(txt, " Temperature (mK)");
  //    fprintf(file, "%6.2f %6.2f moveto\n 90 rotate\n (%s) show\n -90 rotate\n",xoffset -52.0, 265.0, txt);
  fprintf(file, "%6.2f %6.2f moveto\n 90 rotate\n (%s) show\n -90 rotate\n", xoffset - 30.0, 150.0, txt);

  sprintf(txt, "Frequency (MHz)");
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 160.0, 65.0, txt);
  fprintf(file, "showpage\n%c%cTrailer\n", '%', '%');
  fclose(file);
}

/*
double agauss(double frq,double feor,double wid)
{ double a;
if(frq<feor) a= exp(-0.693*(frq-feor)*(frq-feor)/(wid*wid*2.0*0.25));
else a= exp(-0.693*(frq-feor)*(frq-feor)/(wid*wid*0.25));
return a;
}
*/

double agauss(double frq, double feor, double wid) {
  double a, feor1, feor2;
  feor1 = feor + 10;
  feor2 = feor - 10;
  a = 0.5 * (tanh((1420 / frq - 1420 / feor1) / (wid * 1420 / (feor1 * feor1))) + 1.0);
  a += -0.5 * (tanh((1420 / frq - 1420 / feor2) / (wid * 1420 / (feor2 * feor2))) + 1.0);
  return a;
}

void plotq(unsigned char map[], double fcen, double dfcen, double amp, double damp, double wid, double dwid, double tau, double dtau, char rslts[]) {
  FILE *file2;
  double ampp;
  double gray, h, s, b;
  double x1, y1, x2, y2, dx, dy, dd, w, ww;
  double a1, a2, d, max, minn;
  double rlev, sm, mm;
  double mapp[16][16];
  int n, xx, yy, xx1, xx2, k1, k2, yy1, yy2, side;
  int i, j, k, l, m, min;
  int bmode, xoff, yoff, mode;
  //    time_t now;

  for (mode = 1; mode <= 6; mode++) {
    if (mode == 1) {
      for (j = 0; j < 16; j++) {
        for (i = 0; i < 16; i++) {
          min = 99;
          {
            for (k = 0; k < 16; k++)
              for (l = 0; l < 16; l++) {
                m = i + j * 16 + k * 16 * 16 + l * 16 * 16 * 16;
                if (map[m] < min) min = map[m];
              }
          }
          mapp[i][j] = min / 10.0;
        }
      }
      printf("frequency across\n");
      xoff = 100;
      yoff = 600;
    }

    if (mode == 2) {
      for (k = 0; k < 16; k++) {
        for (i = 0; i < 16; i++) {
          min = 99;
          {
            for (j = 0; j < 16; j++)
              for (l = 0; l < 16; l++) {
                m = i + j * 16 + k * 16 * 16 + l * 16 * 16 * 16;
                if (map[m] < min) min = map[m];
              }
          }
          mapp[i][k] = min / 10.0;
        }
      }
      printf("frequency across\n");
      xoff = 300;
      yoff = 600;
    }

    if (mode == 3) {
      for (l = 0; l < 16; l++) {
        for (i = 0; i < 16; i++) {
          min = 99;
          {
            for (j = 0; j < 16; j++)
              for (k = 0; k < 16; k++) {
                m = i + j * 16 + k * 16 * 16 + l * 16 * 16 * 16;
                if (map[m] < min) min = map[m];
              }
          }
          mapp[i][l] = min / 10.0;
        }
      }
      printf("frequency across\n");
      xoff = 100;
      yoff = 400;
    }

    if (mode == 4) {
      for (k = 0; k < 16; k++) {
        for (j = 0; j < 16; j++) {
          min = 99;
          {
            for (i = 0; i < 16; i++)
              for (l = 0; l < 16; l++) {
                m = i + j * 16 + k * 16 * 16 + l * 16 * 16 * 16;
                if (map[m] < min) min = map[m];
              }
          }
          mapp[j][k] = min / 10.0;
        }
      }
      printf("frequency across\n");
      xoff = 300;
      yoff = 400;
    }

    if (mode == 5) {
      for (l = 0; l < 16; l++) {
        for (j = 0; j < 16; j++) {
          min = 99;
          {
            for (i = 0; i < 16; i++)
              for (k = 0; k < 16; k++) {
                m = i + j * 16 + k * 16 * 16 + l * 16 * 16 * 16;
                if (map[m] < min) min = map[m];
              }
          }
          mapp[j][l] = min / 10.0;
        }
      }
      printf("frequency across\n");
      xoff = 100;
      yoff = 200;
    }

    if (mode == 6) {
      for (l = 0; l < 16; l++) {
        for (k = 0; k < 16; k++) {
          min = 99;
          {
            for (i = 0; i < 16; i++)
              for (j = 0; j < 16; j++) {
                m = i + j * 16 + k * 16 * 16 + l * 16 * 16 * 16;
                if (map[m] < min) min = map[m];
              }
          }
          mapp[k][l] = min / 10.0;
        }
      }
      printf("frequency across\n");
      xoff = 300;
      yoff = 200;
    }

    x1 = x2 = y1 = y2 = xx1 = xx2 = yy1 = yy2 = 0;
    bmode = 0;
    if (mode == 1) {
      if ((file2 = fopen("spem.pos", "w")) == NULL) {
        printf(" output error\n");
        return;
      }
      fprintf(file2, "%%!PS-Adobe-3.0 EPSF-3.0\n"); /* %%Trailer - to print this file */
      fprintf(file2, "%%%%BoundingBox:  0 0 612 792\n%%%%EndProlog\n");
      fprintf(file2, "/box1\n{1 0 rlineto\n 0 1 rlineto\n -1 0 rlineto\n closepath } def\n");
    }
    fprintf(file2, "2 setlinewidth\n");
    //    now = time(NULL);
    fprintf(file2, "0 0 0 sethsbcolor\n");
    fprintf(file2, "newpath\n %d %d moveto\n %d %d lineto\nstroke\n", xoff, yoff, xoff + 128, yoff);
    fprintf(file2, "newpath\n %d %d moveto\n %d %d lineto\nstroke\n", xoff, yoff + 128, xoff + 128, yoff + 128);
    fprintf(file2, "newpath\n %d %d moveto\n %d %d lineto\nstroke\n", xoff, yoff, xoff, yoff + 128);
    fprintf(file2, "newpath\n %d %d moveto\n %d %d lineto\nstroke\n", xoff + 128, yoff, xoff + 128, yoff + 128);

    fprintf(file2, "newpath\n %d %d moveto\n %d %d lineto\nstroke\n", xoff + 64, yoff, xoff + 64, yoff - 4);
    fprintf(file2, "newpath\n %d %d moveto\n %d %d lineto\nstroke\n", xoff + 64, yoff + 128, xoff + 64, yoff + 132);
    fprintf(file2, "newpath\n %d %d moveto\n %d %d lineto\nstroke\n", xoff, yoff + 64, xoff - 4, yoff + 64);
    fprintf(file2, "newpath\n %d %d moveto\n %d %d lineto\nstroke\n", xoff + 128, yoff + 64, xoff + 132, yoff + 64);
    fprintf(file2, "/Times-Roman findfont\n 12 scalefont\n setfont\n");
    max = 0;
    minn = 1e99;
    sm = 16;
    mm = 0;
    for (xx = 0; xx < 128; xx += 1) {
      for (yy = 0; yy < 128; yy += 1) {
        w = ww = 0;
        for (i = 0; i < 16; i++) {
          for (j = 0; j < 16; j++) {
            dx = (xx + mm - i * 8);
            dy = (yy + mm - j * 8);
            dd = exp((-dx * dx - dy * dy) / sm);
            w += mapp[i][j] * dd;
            ww += dd;
            //  w +=  mapp[i][j]/dd;
            //  ww += 1.0/dd;
          }
        }
        ampp = w / ww;
        if (ampp > max) max = ampp;
        if (ampp < minn) minn = ampp;
      }
    }
    printf("max %f min %f\n", max, minn);
    rlev = 3;
    if (mode == 1) {
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff + 60, yoff - 14, fcen);
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff + 60 - 64, yoff - 14, fcen - 8 * dfcen);
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff + 60 + 64, yoff - 14, fcen + 8 * dfcen);
      fprintf(file2, "%d %d  moveto\n (freq(MHz)) show\n", xoff + 50, yoff - 24);
      fprintf(file2, "%d %d  moveto\n (%4.2f) show\n", xoff - 30, yoff + 64 + 64, amp + 8 * damp);
      fprintf(file2, "%d %d  moveto\n (%4.2f) show\n", xoff - 30, yoff + 64 - 64, amp - 8 * damp);
      fprintf(file2, "%d %d  moveto\n (%4.2f) show\n", xoff - 30, yoff + 64, amp);
      fprintf(file2, "%d %d  moveto\n 90 rotate\n (amp) show\n -90 rotate\n", xoff - 40, yoff + 64);
    }

    if (mode == 2) {
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff + 60, yoff - 14, fcen);
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff + 60 - 64, yoff - 14, fcen - 8 * dfcen);
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff + 60 + 64, yoff - 14, fcen + 8 * dfcen);
      fprintf(file2, "%d %d  moveto\n (freq(MHz)) show\n", xoff + 50, yoff - 24);
      fprintf(file2, "%d %d  moveto\n (%4.2f) show\n", xoff - 34, yoff + 64 + 64, wid + 8 * dwid);
      fprintf(file2, "%d %d  moveto\n (%4.2f) show\n", xoff - 34, yoff + 64 - 64, wid - 8 * dwid);
      fprintf(file2, "%d %d  moveto\n (%4.2f) show\n", xoff - 34, yoff + 64, wid);
      fprintf(file2, "%d %d  moveto\n 90 rotate\n (width) show\n -90 rotate\n", xoff - 40, yoff + 64);
    }

    if (mode == 3) {
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff + 60, yoff - 14, fcen);
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff + 60 - 64, yoff - 14, fcen - 8 * dfcen);
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff + 60 + 64, yoff - 14, fcen + 8 * dfcen);
      fprintf(file2, "%d %d  moveto\n (freq(MHz)) show\n", xoff + 50, yoff - 24);
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff - 25, yoff + 64 + 64, tau + 8 * dtau);
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff - 25, yoff + 64 - 64, tau - 8 * dtau);
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff - 25, yoff + 64, tau);
      fprintf(file2, "%d %d  moveto\n 90 rotate\n (tau) show\n -90 rotate\n", xoff - 40, yoff + 64);
    }
    if (mode == 4) {
      fprintf(file2, "%d %d  moveto\n (%4.2f) show\n", xoff + 60, yoff - 14, amp);
      fprintf(file2, "%d %d  moveto\n (%4.2f) show\n", xoff + 60 - 64, yoff - 14, amp - 8 * damp);
      fprintf(file2, "%d %d  moveto\n (%4.2f) show\n", xoff + 60 + 64, yoff - 14, amp + 8 * damp);
      fprintf(file2, "%d %d  moveto\n (amp) show\n", xoff + 60, yoff - 24);
      fprintf(file2, "%d %d  moveto\n (%4.2f) show\n", xoff - 34, yoff + 64 + 64, wid + 8 * dwid);
      fprintf(file2, "%d %d  moveto\n (%4.2f) show\n", xoff - 34, yoff + 64 - 64, wid - 8 * dwid);
      fprintf(file2, "%d %d  moveto\n (%4.2f) show\n", xoff - 34, yoff + 64, wid);
      fprintf(file2, "%d %d  moveto\n 90 rotate\n (width) show\n -90 rotate\n", xoff - 40, yoff + 64);
    }

    if (mode == 5) {
      fprintf(file2, "%d %d  moveto\n (%4.2f) show\n", xoff + 60, yoff - 14, amp);
      fprintf(file2, "%d %d  moveto\n (%4.2f) show\n", xoff + 60 - 64, yoff - 14, amp - 8 * damp);
      fprintf(file2, "%d %d  moveto\n (%4.2f) show\n", xoff + 60 + 64, yoff - 14, amp + 8 * damp);
      fprintf(file2, "%d %d  moveto\n (amp) show\n", xoff + 60, yoff - 24);
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff - 25, yoff + 64 + 64, tau + 8 * dtau);
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff - 25, yoff + 64 - 64, tau - 8 * dtau);
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff - 25, yoff + 64, tau);
      fprintf(file2, "%d %d  moveto\n 90 rotate\n (tau) show\n -90 rotate\n", xoff - 40, yoff + 64);
    }
    if (mode == 6) {
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff + 60, yoff - 14, wid);
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff + 60 - 64, yoff - 14, wid - 8 * dwid);
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff + 60 + 64, yoff - 14, wid + 8 * dwid);
      fprintf(file2, "%d %d  moveto\n (width) show\n", xoff + 60, yoff - 24);
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff - 25, yoff + 64 + 64, tau + 8 * dtau);
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff - 25, yoff + 64 - 64, tau - 8 * dtau);
      fprintf(file2, "%d %d  moveto\n (%4.1f) show\n", xoff - 25, yoff + 64, tau);
      fprintf(file2, "%d %d  moveto\n 90 rotate\n (tau) show\n -90 rotate\n", xoff - 40, yoff + 64);
    }

    for (xx = 0; xx < 128; xx += 1) {
      for (yy = 0; yy < 128; yy += 1) {
        w = ww = 0;
        for (i = 0; i < 16; i++) {
          for (j = 0; j < 16; j++) {
            dx = (xx + mm - i * 8);
            dy = (yy + mm - j * 8);
            dd = dx * dx + dy * dy;
            dd = exp((-dx * dx - dy * dy) / sm);
            w += mapp[i][j] * dd;
            ww += dd;
          }
        }
        ampp = (w / ww) / max;
        if (ampp > 1.0) ampp = 1.0;
        h = ampp;
        s = b = 1.0;
        gray = ampp;
        if (ampp > 0.8) {
          h = 0;
          s = 0;
          b = 1;
        }
        h = 0.67 * floor(h * rlev) / rlev;
        if (h < 0.2) h = 0;
        gray = floor(gray * rlev) / rlev;
        if (bmode == 0 && ampp > 0.8 * 0)
          fprintf(file2, "newpath\n %d %d moveto box1\n %3.1f %3.1f %3.1f sethsbcolor fill\n", xx + xoff, yy + yoff, h, s, b);
        if (bmode == 1 && ampp > 0.8) fprintf(file2, "newpath\n %d %d moveto box1\n %3.2f setgray fill\n", xx + xoff, yy + yoff, gray);
      }
    }
    fprintf(file2, "0 0 0 sethsbcolor\n");

    for (xx = 0; xx < 128; xx += 1) {
      for (yy = 0; yy < 128; yy += 1) {
        n = 0;
        for (side = 0; side < 4; side++) {
          if (side == 0) {
            xx1 = xx;
            yy1 = yy;
            xx2 = xx + 1;
            yy2 = yy;
          }
          if (side == 1) {
            xx1 = xx + 1;
            yy1 = yy;
            xx2 = xx + 1;
            yy2 = yy + 1;
          }
          if (side == 2) {
            xx1 = xx;
            yy1 = yy + 1;
            xx2 = xx + 1;
            yy2 = yy + 1;
          }
          if (side == 3) {
            xx1 = xx;
            yy1 = yy;
            xx2 = xx;
            yy2 = yy + 1;
          }
          w = ww = 0;
          for (i = 0; i < 16; i++) {
            for (j = 0; j < 16; j++) {
              dx = (xx2 + mm - i * 8);
              dy = (yy2 + mm - j * 8);
              dd = dx * dx + dy * dy;
              dd = exp((-dx * dx - dy * dy) / sm);
              w += mapp[i][j] * dd;
              ww += dd;
            }
          }
          ampp = (w / ww) / max;
          if (ampp > 1.1) ampp = 1.1;
          a2 = rlev * ampp;
          w = ww = 0;
          for (i = 0; i < 16; i++) {
            for (j = 0; j < 16; j++) {
              dx = (xx1 + mm - i * 8);
              dy = (yy1 + mm - j * 8);
              dd = dx * dx + dy * dy;
              dd = exp((-dx * dx - dy * dy) / sm);
              w += mapp[i][j] * dd;
              ww += dd;
            }
          }
          ampp = (w / ww) / max;
          if (ampp > 1.1) ampp = 1.1;
          a1 = rlev * ampp;
          k1 = a1;
          k2 = a2;
          if (k2 != k1) {
            n++;
            if (a2 > a1)
              d = (k2 - a1) / (a2 - a1);
            else
              d = (a1 - k1) / (a1 - a2);
            if (side == 0) {
              x1 = xx1 + d;
              y1 = yy1;
            }
            if (side == 1 && n == 1) {
              x1 = xx1;
              y1 = yy1 + d;
            }
            if (side == 1 && n == 2) {
              x2 = xx1;
              y2 = yy1 + d;
            }
            if (side == 2 && n == 1) {
              x1 = xx1 + d;
              y1 = yy1;
            }
            if (side == 2 && n == 2) {
              x2 = xx1 + d;
              y2 = yy1;
            }
            if (side == 3) {
              x2 = xx1;
              y2 = yy1 + d;
            }
          }
        }
        if (n == 2) fprintf(file2, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\nstroke\n", x1 + xoff, y1 + yoff, x2 + xoff, y2 + yoff);
      }
    }
    fprintf(file2, "20 50 moveto\n (%s) show\n", rslts);
  }
  fprintf(file2, "showpage\n%%%%Trailer\n");
  fclose(file2);
  return;
}
