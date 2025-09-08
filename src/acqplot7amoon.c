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
#define NFIT 200
double polyfit(int, int, double *, double *, double *, double *, int);
double fitfun(int, int, int);
double rmscalc(int, double *, double *);
void rmsfilt(int, int, double, double *, double *);
double gst(double);
void sunradec(double, double *, double *);
void moonradec(double, double *, double *, double, double);
void radec_azel(double, double, double, double *, double *);
double tosecs(int, int, int, int, int);
void toyrday(double, int *, int *, int *, int *, int *);
void dsmooth(int, int, double *, double *);
double tempmon(double, double, double *);
void plotfspec(int, double, double *, double *, int);
void waterspec(int, int, double *, int, int, int, int, double, double);
void waterwrite(int, int, double *, int, int, int, int, int, double);
void specout(int, double *, double *, char *, int, int, int, int, int, int, double, double, double, double, double);
void qrd(long double *, int, double *);
void MatrixInvert(int);
void b64init(int *);
void rmsinit(double *);
double rmsf(double *);
static int mode = 0;
static int pfit = 0;
static double scaledb = 0;
// static long double aarr[100000];
// static double bbrr[10000];
// static long double arrr[100000];
static double mcalc[65536 * 40];
static double data[100000];
static double data2[100000];
static double avdata[100000];
static double avdata0[100000];
static double avdata1[100000];
static double avdata2[100000];
static double specm[100000], specr[100000], specmav[100000], specrav[100000], wttav[100000], spec[100000], specav[100000];
static double avd0[100000], avd1[100000], avd2[100000];
static double sploadp[100000], spload[100000];
static double wtt[100000];
static double wtts[100000];
static double wtts2[100000];
static double avdataw[100000];
static double dataout[100000];

static double ncdata[100000];
static double nndata[100000];
static double ffun[100000];

char fname[256];
int tstart, tstop, smooth, nline;
double fstart0, fstart, fstep, fstop, fstop0, secint;
double freqstart, freqstop, rfi, resol, peakpwr, ndataw, pkpwrm, wttt, maxscale, minscale;
int main(int argc, char *argv[]) {
  static char buf[1000000];
  int i, j, k, nrfi, m, line, np, nnp, day, pp0, pp1, wwd, zwt, zwtp, nspec, temp, nblock, jj, ns, nns, nt, ntt, nn;
  int yr, dy, hr, mn, sc, swpos, trec, nstart, water, ppmax, recom, fmode;
  int b64[256];
  unsigned int kk, v24;  // for speedy sscanf
  double avpwr, avpwr2, avpwrt, pwr, tcal, tload, p0, p1, p2, rms, rmss, secs, delaystart, secsst;
  double f, rmstheory, pkpwr, adcov, padcov, gal, gha, dgha, ghav, ghamax, ghamin, avsec, orbpwr, orbpwr2, sweep, maxsweep, mxsweep;
  double avpkpwr, a, b, ppercent, lst, ha, ra, dec, az, el, lat, lon, sunlim, moonlim, tav, tr, tempr, wscalemax, wscalemin, max, min, av, d150,
      dload, dloadmax, minpwr, maxrmsf, rrmsf;
  double a00, a01, a10, a11, b0, b1, d, aa00, aa10, aa01, aa11, n0, n1, n2, maxp0, minp0, avp0, np0, mel, maxfm, fmpwr;
  FILE *file1;
  FILE *timeflag_file;
  char *p;
  short day_flag, gal_flag, gha_flag, nstart_flag, sunlim_flag, moonlim_flag, delaystart_flag, 
      adc_flag, buf_flag, tstop_hr_flag, ppercent_flag, pkpower_flag, d150_flag, 
      dloadmax_flag, rmsf_flag, fmpwr_flag;
  b64init(b64);
  rmsinit(ffun);
  fstart = f = 0.0;
  tstart = 0;
  tstop = 24;
  day = 0;
  trec = 0;
  wwd = 0;
  recom = 0;
  fmode = 0;
  fstep = 0.15625;
  tcal = 400;
  tload = 300;
  tav = tload;
  tr = 60;  // trec
  rfi = 0;
  nrfi = 0;
  nline = 0;
  ns = -1;
  nns = nt = ntt = 0;
  nn = 256;
  d150 = 1e99;  // 20 best
  smooth = 0;
  water = 0;
  sunlim = 1e99;
  moonlim = 0;
  padcov = 1;
  avsec = 0;
  ppmax = 0;
  minpwr = 0;
  maxp0 = -1e99;
  minp0 = 1e99;
  np0 = avp0 = 1e-99;
  gal = 0;
  gha = 0;
  dgha = 48;
  peakpwr = 1e99;
  pkpwrm = 1e99;
  maxsweep = 0;
  mxsweep = 1e99;
  maxrmsf = 1e99;
  freqstart = 40.0;
  freqstop = 200.0;
  nspec = 16384;
  nstart = -1;
  maxscale = minscale = 0;
  wscalemax = 1000;
  wscalemin = 0;
  secsst = 0;
  delaystart = -1e99;
  dloadmax = 1e99;
  maxfm = 1e99;
  dload = 0;
  m = 0;
  if (argc > 1) sscanf(argv[1], "%s", fname);
  for (i = 0; i < argc; i++) {
    sscanf(argv[i], "%79s", buf);
    if (strstr(buf, "-mode")) { sscanf(argv[i + 1], "%d", &mode); }
    if (strstr(buf, "-fmode")) { sscanf(argv[i + 1], "%d", &fmode); }
    if (strstr(buf, "-scaledb")) { sscanf(argv[i + 1], "%lf", &scaledb); }
    if (strstr(buf, "-wscalemax")) { sscanf(argv[i + 1], "%lf", &wscalemax); }
    if (strstr(buf, "-wscalemin")) { sscanf(argv[i + 1], "%lf", &wscalemin); }
    if (strstr(buf, "-tstart")) { sscanf(argv[i + 1], "%d", &tstart); }
    if (strstr(buf, "-tstop")) { sscanf(argv[i + 1], "%d", &tstop); }
    if (strstr(buf, "-pfit")) { sscanf(argv[i + 1], "%d", &pfit); }
    if (strstr(buf, "-ns")) { sscanf(argv[i + 1], "%d", &ns); }
    if (strstr(buf, "-nn")) { sscanf(argv[i + 1], "%d", &nn); }
    if (strstr(buf, "-d150")) { sscanf(argv[i + 1], "%lf", &d150); }
    if (strstr(buf, "-dloadmax")) { sscanf(argv[i + 1], "%lf", &dloadmax); }
    if (strstr(buf, "-smooth")) { sscanf(argv[i + 1], "%d", &smooth); }
    if (strstr(buf, "-day")) { sscanf(argv[i + 1], "%d", &day); }
    if (strstr(buf, "-gal")) { sscanf(argv[i + 1], "%lf", &gal); }
    if (strstr(buf, "-gha")) { sscanf(argv[i + 1], "%lf", &gha); }
    if (strstr(buf, "-dgha")) { sscanf(argv[i + 1], "%lf", &dgha); }
    if (strstr(buf, "-nstart")) { sscanf(argv[i + 1], "%d", &nstart); }
    if (strstr(buf, "-tcal")) { sscanf(argv[i + 1], "%lf", &tcal); }
    if (strstr(buf, "-fstart")) { sscanf(argv[i + 1], "%lf", &freqstart); }
    if (strstr(buf, "-fstop")) { sscanf(argv[i + 1], "%lf", &freqstop); }
    if (strstr(buf, "-peakpwr")) { sscanf(argv[i + 1], "%lf", &peakpwr); }
    if (strstr(buf, "-minpwr")) { sscanf(argv[i + 1], "%lf", &minpwr); }
    if (strstr(buf, "-pkpwrm")) { sscanf(argv[i + 1], "%lf", &pkpwrm); }
    if (strstr(buf, "-sunlim")) { sscanf(argv[i + 1], "%lf", &sunlim); }
    if (strstr(buf, "-moonlim")) { sscanf(argv[i + 1], "%lf", &moonlim); }
    if (strstr(buf, "-trec")) { trec = 1; }
    if (strstr(buf, "-water")) { sscanf(argv[i + 1], "%d", &water); }
    if (strstr(buf, "-nrfi")) { sscanf(argv[i + 1], "%d", &nrfi); }
    if (strstr(buf, "-recom")) { sscanf(argv[i + 1], "%d", &recom); }
    if (strstr(buf, "-adcov")) { sscanf(argv[i + 1], "%lf", &padcov); }
    if (strstr(buf, "-maxscale")) { sscanf(argv[i + 1], "%lf", &maxscale); }
    if (strstr(buf, "-minscale")) { sscanf(argv[i + 1], "%lf", &minscale); }
    if (strstr(buf, "-mxsweep")) { sscanf(argv[i + 1], "%lf", &mxsweep); }
    if (strstr(buf, "-maxrmsf")) { sscanf(argv[i + 1], "%lf", &maxrmsf); }
    if (strstr(buf, "-maxfm")) { sscanf(argv[i + 1], "%lf", &maxfm); }
    if (strstr(buf, "-delaystart")) { sscanf(argv[i + 1], "%lf", &delaystart); }
    if (strstr(buf, "-wwd")) { wwd = 1; }
    if (strstr(buf, "-rfi")) { sscanf(argv[i + 1], "%lf", &rfi); }
  }
  //  if(water && smooth > 4) water = 2;
  for (i = 0; i < 100000; i++) {
    data[i] = avdata[i] = avdata0[i] = avdata1[i] = avdata2[i] = ncdata[i] = nndata[i] = 0;
    spload[i] = sploadp[i] = 0;
    spec[m] = specav[m] = specm[i] = specr[i] = specmav[i] = specrav[i] = wttav[i] = 0;
    avd0[i] = avd1[i] = avd2[i] = avdataw[i] = 0;
    wtt[i] = 1.0;
    wtts[i] = 1.0;
  }
  // corrected formula from memo 174
  n0 = tav + tr;
  n1 = -tr - tload + (tav - tload) * (tload + tr) / tcal;
  n2 = -tav + tload - (tav - tload) * (tload + tr) / tcal;
  rmstheory = sqrt(n0 * n0 + n1 * n1 + n2 * n2);
  pp0 = pp1 = 0;
  lst = 0;
  np = 0;
  avpwr = 0;
  rms = rmss = 0;
  ndataw = 0;
  ghav = 0;
  ghamax = -120;
  ghamin = 120;
  ppercent = 0;
  avpkpwr = 0;
  zwt = 0;
  nnp = 1;
  j = 0;
  if ((file1 = fopen(fname, "r")) == NULL) {
    printf("cannot open file:%s\n", fname);
    return 0;
  }
  line = nspec = 0;
  printf("file %s\n", fname);

  if ((timeflag_file = fopen("time_flags.txt", "w")) == NULL) {
    printf("cannot open file: time_flags.txt\n");
    return 0;
  }

  fprintf(timeflag_file, 
    "# HA0 day0 gal0 gha0 nstart0 sunlim0 moonlim0 delaystart0 adc0 buf0 tstop-hr0 "
    "HA1 day1 gal1 gha1 nstart1 sunlim1 moonlim1 delaystart1 adc1 buf1 tstop-hr1 "
    "HA2 day2 gal2 gha2 nstart2 sunlim2 moonlim2 delaystart2 adc2 buf2 tstop-hr2 "
    "ppercent pkpower d150 dloadmax rmsf fmpwr\n"
  );
  long foffset;
  int prev_swpos = -1;

  foffset = ftell(timeflag_file);
  
  while (fgets(buf, 1000000, file1) != 0) {
    // Record the start of each line in the timeflag file.
    

    if (buf[0] == '#' && buf[2] == 's') sscanf(buf, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %d %*s %*s %d %*s %d", &temp, &nblock, &nspec);  // db310
    // skip this line if it is a comment
    if (buf[0] == '*' || buf[0] == '\n' || buf[0] == '#' || buf[0] == ';') continue;
    
    line++;
    avpwr = avpwr2 = avpwrt = 0;
    sscanf(buf, "%4d:%3d:%2d:%2d:%2d %d %lf %lf %lf %lf", &yr, &dy, &hr, &mn, &sc, &swpos, &fstart, &fstep, &fstop, &adcov);
    if (freqstart < fstart) freqstart = fstart;
    if (fabs(fstep - 0.006104) < 0.001) fstep = 200.0 / nspec;  // for higher accuracy
    if (nspec == 0) {
      if (fstep > 0.006408 && fstep < 0.006410) fstep = 420.0 / 65536.0;
      if (fstep > 0.012)
        nspec = 16384;
      else {
        if (fstep > 0.006)
          nspec = 32768;
        else
          nspec = 65536;
      }
    }
    lon = 116.5 * PI / 180.0;
    lat = -26.7 * PI / 180.0;
    secs = tosecs(yr, dy, hr, mn, sc);
    if (!secsst) secsst = secs;
    ha = (gst(secs) - (17.0 + 45.67 / 60.0) * PI / 12.0 + lon) * 12.0 / PI;
    if (ha > 12.0) ha += -24.0;
    if (ha <= -12.0) ha += 24.0;
    if (sunlim < 1e3) {
      sunradec(secs, &ra, &dec);
      radec_azel(gst(secs) - ra + lon, dec, lat, &az, &el);
      el = el * 180.0 / PI;
    } else
      el = 0;
    if (moonlim) {
      moonradec(secs, &ra, &dec, lat, lon);
      radec_azel(gst(secs) - ra + lon, dec, lat, &az, &mel);
      mel = mel * 180.0 / PI;
    } else
      mel = 0;
    lst = (gst(secs) + lon) * 12.0 / PI;
    if (lst > 24.0) lst += -24.0;
    if (lst < 0.0) lst += 24.0;
    printf("%4d:%03d:%02d:%02d:%02d %5.2f %1.10e %.2f %1.10e ", yr, dy, hr, mn, sc, lst, ha, secs, gst(secs));
    j = 0;

    day_flag = (dy == day || day == 0);
    gal_flag =(gal == 0 || (gal > 0.0 && fabs(ha) < gal) || (gal < 0.0 && fabs(ha) >= 6.0));
    gha_flag =(fabs(gha - ha) < dgha || fabs(gha - ha - 24.0) < dgha || gha > 1e3);
    nstart_flag = (nstart == -1 || (nstart <= line && nstart > 0)) && (nstart != -999);
    sunlim_flag =el < sunlim;
    moonlim_flag=((mel < moonlim && moonlim < 0) || (mel > moonlim && moonlim > 0) || !moonlim);
    delaystart_flag=(secs - secsst > delaystart); 
    adc_flag=adcov < padcov;
    buf_flag=strlen(buf) > 131000;
    tstop_hr_flag=((hr >= tstart && hr <= tstop) || (tstart > tstop && (hr <= tstop || hr >= tstart)));

    // Before we start writing to the timeflag file, check if this swpos is correct.
    if ((swpos != (prev_swpos + 1)) && !(swpos==0 && prev_swpos == 2)){
      fseek(timeflag_file, foffset, SEEK_SET);
    }
    prev_swpos = swpos;

    fprintf(timeflag_file, "%lf %d %d %d %d %d %d %d %d %d %d ", 
      ha, day_flag, gal_flag, gha_flag, nstart_flag, sunlim_flag, moonlim_flag, 
      delaystart_flag, adc_flag, buf_flag, tstop_hr_flag
    );

    if (day_flag && gal_flag && gha_flag && nstart_flag && sunlim_flag && moonlim_flag && delaystart_flag && adc_flag && buf_flag && tstop_hr_flag)
    {
      // Everything above got through fine.
      p = buf;
      p = strstr(p, "spectrum");
      if (*p && p) p = strchr(p, ' ');
      while (*p == ' ') p++;
      j = 0;
      pkpwr = 0;
      tempr = 0;
      orbpwr = orbpwr2 = 0;
      while (j < nspec && p && *p != '\n') {
        kk = b64[(int)*p];
        p++;
        v24 = (kk << 18);
        kk = b64[(int)*p];
        p++;
        v24 |= (kk << 12);
        kk = b64[(int)*p];
        p++;
        v24 |= (kk << 6);
        kk = b64[(int)*p];
        p++;
        v24 |= kk;
        pwr = -((double)v24) * 1e-5;
        if (swpos == 0) {
          avdata0[j] = pwr;
          pp0 = line;
        }
        if (swpos == 1) {
          avdata1[j] = pwr;
          pp1 = line;
        }
        if (swpos == 2 && (pp0 == line - 2) && (pp1 == line - 1)) {  // need to get adjacent swpos
          avdata2[j] = pwr;
          p0 = pow(10.0, 0.1 * avdata0[j]);
          p1 = pow(10.0, 0.1 * avdata1[j]);
          p2 = pow(10.0, 0.1 * avdata2[j]);
          avd0[j] += p0;
          avd1[j] += p1;
          avd2[j] += p2;
          if (tcal > 0) {
            if (p2 - p1 > 1e-99)
              data[j] = ((p0 - p1) / (p2 - p1)) * tcal + tload;
            else
              data[j] = 0;
            spload[j] = ((p1) / (p2 - p1)) * tcal + tload;
          } else
            data[j] = (-(p0 - p1) / (p2)) * tcal + tload;
          if (trec) data[j] = ((p1) / (p2 - p1)) * tcal - tload;
          if (mode == 1) data[j] = p0 * 300.0e08;
          if (mode == 2) data[j] = p1 * 300.0e08;
          if (mode == 3) data[j] = p2 * 300.0e08;
          if (mode == 4) data[j] = (p2 / p1) * tcal;
          if (mode == 5) data[j] = tcal * (p2 - p1) / p1;
          if (mode == 6) data[j] = tcal * (p0 - p2) / (p2 - p1) + tcal + tload;
          if (mode == 7) data[j] = tcal * (p1 - p0) / p0;
          if (!mode && j < (1000 * nspec) / 32768) data[j] = 0;
          //                    printf("j %d v24 %d pwr %f data %f\n",j,v24,pwr,data[j]);
          if (data[j] > pkpwr && fstart + j * fstep > 80.0) pkpwr = data[j];
          if (data[j] > orbpwr && fabs(fstart + j * fstep - 137.5) < 1.0) orbpwr = data[j];
        }
        avpwr += pwr;
        if (fstart + j * fstep > 100.0) avpwr2 += 50.0 * pow(10.0, 0.1 * pwr);
        avpwrt += 50.0 * pow(10.0, 0.1 * pwr);
        j++;
      }
      if (swpos == 0) { ppercent = 100.0 * avpwr2 / avpwrt; }
      if (swpos == 2 && (pp0 == line - 2) && (pp1 == line - 1)) {
        pp0 = pp1 = 0;

        maxsweep = 1e-99;
        for (i = 2; i < j - 2; i++) {
          if (fstart + i * fstep >= 50.0) {
            a = b = 0;
            for (m = 0; m < 1; m++) {
              a += data[i + 1 + m];
              b += data[i - 1 - m];
            }
            a = (a + b) / 2.0;
            sweep = (data[i] - a);
            if (fabs(sweep) > maxsweep) maxsweep = fabs(sweep);
          }
        }

        rrmsf = rmsf(data);

        fmpwr = 0;
        for (i = 0; i < j; i++) {
          if (fstart + i * fstep > 88 && fstart + i * fstep < 120) {
            a = data[i] - 0.5 * (data[i + 1] + data[i - 1]);
            if (a > fmpwr) fmpwr = a;
          }
        }

        a = b = 1e-99;
        for (i = 0; i < j; i++) {
          if (fabs(fstart + i * fstep - 137.9) < 0.02) orbpwr2 = data[i] - data[i + 50];
          if (data[i] > 0.0 && data[i] < pkpwr / 10.0 && fstart + i * fstep > 80.0) {
            a += data[i];
            b++;
          }
        }
        pkpwr = 10.0 * log10(pkpwr / (a / b));
        orbpwr = 10.0 * log10(orbpwr / (a / b));
        if (orbpwr > pkpwr) pkpwr = orbpwr;
        pkpwr = orbpwr;  // O.K. orbpwr and peak power sufficient for verystrong signals
        if (orbpwr2 > 0) orbpwr2 = 10.0 * log10(orbpwr2 / (a / b));
        //                if(orbpwr2 > 5.0) pkpwr = 10.0*orbpwr2;
        tempr = tempmon(fstart, fstep, avdata0);

        if (dloadmax < 1e6) {
          av = 0;
          a = 0;
          for (i = 0; i < j; i++) {
            if (fstart + i * fstep > 50.0 && fstart + i * fstep < 198.0) {
              av += sploadp[i] - spload[i];
              a++;
            }
          }
          if (sploadp[j / 2])
            dload = fabs(av / a) * sqrt(nblock * (double)nspec * 100e6 / 400e6) / rmstheory;  // fix for fastspec
          else
            dload = 0;
          if (dload > dloadmax) printf("dload %f\n", dload);
          for (i = 0; i < j; i++) sploadp[i] = spload[i];
        }

        d = 0;
        if (d150 < 1e6) {
          av = a = 0.0;
          for (i = 0; i < j; i++) {
            if (fabs(fstart + i * fstep - 153.5) < 1.5) {
              av += data[i];
              a++;
            }
          }
          rms = 0;
          av = av / a;
          for (i = 0; i < j; i++) {
            if (fabs(fstart + i * fstep - 153.5) < 1.5) { rms += (data[i] - av) * (data[i] - av); }
          }
          av = b = 0;
          for (i = 0; i < j; i++) {
            if (fabs(fstart + i * fstep - 157.0) < 1.5) {
              av += data[i];
              b++;
            }
          }
          d = 200.0 * sqrt(rms / a) / (av / b);
          if (d >= d150) printf("d150 %f %f\n", d, av / b);
        }

        if ((ppercent < peakpwr && (ppercent > minpwr)) && ((pkpwr < pkpwrm) || (pkpwr > fabs(pkpwrm) && pkpwrm < 0)) && d <= d150 &&
            dload < dloadmax && rrmsf < maxrmsf && fmpwr < maxfm && j == nspec)
        {
          fprintf(timeflag_file, " 1 1 1 1 1 1\n");
          foffset = ftell(timeflag_file);
  
          wttt = 1;
          if (ppercent > maxp0) maxp0 = ppercent;
          if (ppercent < minp0) minp0 = ppercent;
          avp0 += ppercent;
          np0++;

          for (i = 2; i < j - 2; i++) {
            if (fstart + i * fstep >= 50.0) {
              a = b = 0;
              for (m = 0; m < 1; m++) {
                a += data[i + 1 + m];
                b += data[i - 1 - m];
              }
              a = (a + b) / 2.0;
              sweep = (data[i] - a);
              if (sweep > mxsweep) wtts[i] = 0;
            }
          }

          if (nns < ns) {
            for (i = 0; i < j; i++) data2[i] += data[i];
            nns++;
            nt++;
            ntt++;
          }

          if (nns == ns && nns) {
            zwt = 1;
            zwtp = 0;
            for (i = 0; i < j; i++) wtt[i] = 1;
            for (jj = 0; jj < 10 && zwt > zwtp; jj++) {
              for (i = 0; i < j; i++) {  // sliding rms is better
                if (fstart + i * fstep >= freqstart && fstart + i * fstep <= freqstop) {
                  a00 = a11 = a01 = b0 = b1 = 0;
                  for (k = i - nn; k <= i + nn; k++) {
                    if (k >= 0 && k < j && abs(k - i) > nn / 2) {
                      f = fstart + k * fstep;
                      a00 += wtt[k];
                      a11 += wtt[k] * f * f;
                      a01 += wtt[k] * f;
                      b0 += wtt[k] * data2[k];
                      b1 += wtt[k] * data2[k] * f;
                    }
                  }
                  a10 = a01;
                  d = a00 * a11 - a10 * a01;
                  aa00 = a11 / d;
                  aa10 = -a01 / d;
                  aa01 = -a10 / d;
                  aa11 = a00 / d;
                  a = aa00 * b0 + aa01 * b1;  // constant
                  b = aa10 * b0 + aa11 * b1;  // slope
                  a00 = 0;
                  a11 = 1e-99;
                  for (k = i - nn; k <= i + nn; k++) {
                    if (k >= 0 && k < j && abs(k - i) > nn / 2) {
                      f = fstart + k * fstep;
                      av = a + b * f;
                      a00 += (data2[k] - av) * (data2[k] - av) * wtt[k];
                      a11 += wtt[k];
                    }
                  }
                  if (d) {
                    f = fstart + i * fstep;
                    av = a + b * f;
                    rms = sqrt(a00 / a11);
                    b = (data2[i] - av) / (fabs(rfi * 1.2) * rms);  // for rogue ns = 20 rfi = 2.5 scale up to 3.0
                    if (b > 1)
                      for (k = i - 1; k <= i + 1; k++) wtt[k] = 0;
                    specm[i] = av;
                    spec[i] = (data2[i] - av);
                    specr[i] = wtt[i] * (data2[i] - av);
                  }
                }
              }
              zwtp = zwt;
              m = 0;
              for (i = 0; i < j; i++) {
                if (wtt[i] == 0 && fstart + i * fstep > freqstart && fstart + i * fstep < freqstop) { m++; }
              }
              zwt = m;
              printf("%d zwt %d\n", j, zwt);
            }
            for (i = 0; i < j; i++) {
              data2[i] = 0;
              wttav[i] += nt * wtt[i];
              specmav[i] += specm[i];
              specav[i] += spec[i];
              specrav[i] += specr[i];
            }
            nns = nt = 0;
          }

          for (i = 0; i < j; i++) {
            avdata[i] += wttt * data[i];
            nndata[i] += wttt;
          }
          if (nstart > 0) nstart = -999;
          if (water && nline == 0) {
            waterspec(0, np, avdataw, dy, hr, mn, sc, wscalemax, wscalemin);
            waterwrite(0, np, avdataw, yr, dy, hr, mn, sc, lst);
          }
          for (i = 0; i < nspec; i++) avdataw[i] += data[i];
          avpkpwr += ppercent;
          avsec += secs;
          ndataw++;
          a = ha;
          if (dgha == 48) {
            gha = a;
            dgha = 24;
          }  // needed to get correct ghamax and ghamin when using tstart tstop
          if (a - gha > 12.0) a += -24.0;
          if (a - gha < -12.0) a += 24.0;
          ghav += a;
          a = ha - gha;
          if (a > 12.0) a += -24.0;
          if (a <= -12.0) a += 24.0;
          if (gha + a > ghamax) ghamax = gha + a;
          if (gha + a < ghamin) ghamin = gha + a;

          if (water && (nline % (water) == (water - 1) || peakpwr < 0)) {
            waterwrite(1, np, avdataw, yr, dy, hr, mn, sc, lst);
            waterspec(1, np, avdataw, dy, hr, mn, sc, wscalemax, wscalemin);
          }
          nline++;
        } // if power/rms filters passed
        else{
          fprintf(timeflag_file, "%d %d %d %d %d %d\n", 
            (ppercent < peakpwr && (ppercent > minpwr)), 
            ((pkpwr < pkpwrm) || (pkpwr > fabs(pkpwrm) && pkpwrm < 0)),
            d <= d150,
            dload < dloadmax,
            rrmsf < maxrmsf,
            fmpwr < maxfm
          );
          foffset = ftell(timeflag_file);
        }
      } // if swpos=2 and other swpos got through aux
      else if (swpos == 2) {
        fprintf(timeflag_file, "-1 -1 -1 -1 -1 -1\n");
        foffset = ftell(timeflag_file);
      }
      avpwr = avpwr / (j + 1e-99);
      if (rrmsf < maxrmsf)
        jj = 0;
      else
        jj = 1;
      //        printf("line %4d t %d tot_pwr %6.3f sig_percnt %5.2f cyc_accepted %4d peak %3.0f dB ppmax %d adcov %3.1f swpos %d wttt %2.0f temp
      //        %5.1f orbpwr2 %5.1f rmsf%d %5.1f\n",
      //                       line, temp, avpwrt,100.0*avpwr2/avpwrt,nline,pkpwr,ppmax,adcov,swpos,wttt,tempr,orbpwr2,jj,rrmsf);
      printf(
          "line %4d t %d tot_pwr %6.3f sig_percnt %5.2f cyc_accepted %4d peak %3.0f dB ppmax %d adcov %4.2f swpos %d wttt %2.0f temp %5.1f fmpwr "
          "%5.1f rmsf%d %5.1f\n",
          line, temp, avpwrt, 100.0 * avpwr2 / avpwrt, nline, pkpwr, ppmax, adcov, swpos, wttt, tempr, fmpwr, jj, rrmsf);
      np = j;

    } // endif aux filters passed
    else
    {
      printf(
        "line %4d t %d tot_pwr %6.3f sig_percnt %2d cyc_accepted %4d peak %3.0f dB "
        "ppmax %d adcov %4.2f swpos %d wttt %2d temp %5.1f fmpwr "
        "%5.1f rmsf%d %5.1f",
        line, temp, avpwrt, 0, nline, pkpwr, ppmax, adcov, swpos, 0, tempr, fmpwr, jj, rrmsf
      );
      if(swpos == 2){
        fprintf(timeflag_file, "-1 -1 -1 -1 -1 -1\n");
        foffset = ftell(timeflag_file);
      }
    }
    if (j == 0) printf("\n");
  
  } // while loop over lines in ACQ file
  fclose(file1);

  // If we had a trailing swpos
  printf("ENDED WITH swpos=%d\n", swpos);
  if (swpos != 2){
    fseek(timeflag_file, foffset, SEEK_SET);
    printf("YES I AM TRUNCATING! %ld\n", foffset);
    ftruncate(fileno(timeflag_file), foffset);
  }
  fclose(timeflag_file);

  if (!np) return 0;
  if (nt) { // this is not done for the 2018 Nature paper (i.e. nt=0)
    zwt = 1;
    zwtp = 0;
    for (i = 0; i < np; i++) wtt[i] = 1;
    for (jj = 0; jj < 10 && zwt > zwtp; jj++) {
      for (i = 0; i < np; i++) {  // sliding rms is better
        if (fstart + i * fstep >= freqstart && fstart + i * fstep <= freqstop) {
          a00 = a11 = a01 = b0 = b1 = 0;
          for (k = i - nn; k <= i + nn; k++) {
            if (k >= 0 && k < np && abs(k - i) > nn / 2) {
              f = fstart + k * fstep;
              a00 += wtt[k];
              a11 += wtt[k] * f * f;
              a01 += wtt[k] * f;
              b0 += wtt[k] * data2[k];
              b1 += wtt[k] * data2[k] * f;
            }
          }
          a10 = a01;
          d = a00 * a11 - a10 * a01;
          aa00 = a11 / d;
          aa10 = -a01 / d;
          aa01 = -a10 / d;
          aa11 = a00 / d;
          a = aa00 * b0 + aa01 * b1;  // constant
          b = aa10 * b0 + aa11 * b1;  // slope
          a00 = 0;
          a11 = 1e-99;
          for (k = i - nn; k <= i + nn; k++) {
            if (k >= 0 && k < np && abs(k - i) > nn / 2) {
              f = fstart + k * fstep;
              av = a + b * f;
              a00 += (data2[k] - av) * (data2[k] - av) * wtt[k];
              a11 += wtt[k];
            }
          }
          if (d) {
            f = fstart + i * fstep;
            av = a + b * f;
            rms = sqrt(a00 / a11);
            b = (data2[i] - av) / (fabs(rfi * 1.2) * rms);  // for rogue ns = 20 rfi = 2.5 scale up to 3.0
            if (b > 1)
              for (k = i - 1; k <= i + 1; k++) wtt[k] = 0;
            specm[i] = av;
            spec[i] = (data2[i] - av);
            specr[i] = wtt[i] * (data2[i] - av);
          }
        }
      }
      zwtp = zwt;
      m = 0;
      for (i = 0; i < np; i++) {
        if (wtt[i] == 0 && fstart + i * fstep > freqstart && fstart + i * fstep < freqstop) { m++; }
      }
      zwt = m;
      printf("%d zwt %d\n", np, zwt);
    }
    for (i = 0; i < np; i++) {
      data2[i] = 0;
      wttav[i] += nt * wtt[i];
      specmav[i] += specm[i];
      specav[i] += spec[i];
      specrav[i] += specr[i];
    }
    nns = nt = 0;
  }

  file1 = fopen("avg_temp_and_nsamples.txt", "w");
  for (m = 0; m < np; m++) {
    data[m] = avdata[m] / nndata[m];
    fprintf(file1, "%f %f %f %f\n", fstart + m * fstep, avdata[m], nndata[m], data[m]);

    if (ns > 0) { // not done in B18 (i.e. ns=0)
      data[m] = specmav[m] / ntt;
      if (wttav[m]) data[m] += specrav[m] / wttav[m];
      //          if(wttav[m] == ntt) wtts[m] = 1; else{
      if (wttav[m] == ntt && fabs(fstart + m * fstep - 150.74) > 0.10)
        wtts[m] = 1;
      else {  // weak RFI from geostationary sat
        wtts[m] = 0;
        data[m] = (specmav[m] + specav[m]) / ntt;  // needed to show downweighted points
      }
    }
    if (scaledb > 0) data[m] = 10.0 * log10(fabs(data[m]));
    if (wwd) data[m] = ((avd0[m] - avd1[m]) / (avd2[m] - avd1[m])) * tcal + tload;
    if (recom) data[m] = ((avd0[m]) / avd0[np / 2]) * tcal + tload;
  }
  fclose(file1);

  j = k = 0;
  fstart0 = freqstart;
  fstop0 = freqstop;
  for (i = 0; i < np; i++) { // this is just setting the used frequency data to the front of the array
    f = fstart + i * fstep;
    if (f >= fstart0 && f <= fstop0) {
      if (k == 0) k = i;
      data2[j] = data[i];
      wtts2[j] = wtts[i];
      wtt[j] = 1.0;
      j++;
    }
  }
  np = j;
  //             fstart = fstart0;
  fstart = k * fstep;
  if (freqstop > fstart + (np - 1) * fstep) freqstop = fstart + (np - 1) * fstep;
  for (i = 0; i < np; i++) data[i] = data2[i];
  for (i = 0; i < np; i++) wtts[i] = wtts2[i];
  for (i = 0; i < np; i++) {
    f = fstart + i * fstep;
    if (f >= freqstart && f <= freqstop)
      wtt[i] = 1.0;
    else
      wtt[i] = 0;
  }

  if (pfit) { // this IS done for B18 (pfit=37)
    if (fabs(rfi) > 0.0) { // rfi = 2.5 for B18
      file1 = fopen("rfi_models.txt", "w");
      FILE *file2 = fopen("rfi_flags.txt", "w");
      FILE *file5 = fopen("rfi_windowed_avg.txt", "w");
      FILE *file3 = fopen("rfi_rmss.txt", "w");
      FILE *file4 = fopen("rfi_initial_weights.txt", "w");
      FILE *file6 = fopen("rfi_residuals.txt", "w");
      FILE *file7 = fopen("rfi_sumweights.txt", "w");
      

      zwt = 1;
      zwtp = 0;
      for (i = 0; i < np; i++)
        if (wtts[i] == 0) wtt[i] = 0;
      
      for (j = 0; j < 100 && zwt > zwtp; j++) {
        tav = polyfit(pfit, np, data, mcalc, wtt, dataout, fmode);
        rmss = rmscalc(np, dataout, wtt);

        for(i=0;i<np;i++){
          fprintf(file1, "%1.15e ", data[i] - dataout[i]);
          fprintf(file6, "%1.15e ", dataout[i]);
        }
        fprintf(file1, "\n");
        fprintf(file6, "\n");

        float rmsmax = 0;
        float rmsmin = 0;
        float resmax = 0;
        float resmin = 0;
        float zmax = 0;
        float zmin = 0;

        if (j==0){
          for(i=0; i<np; i++){
            fprintf(file4, "%f\n", wtt[i]);
          }
          fclose(file4);
        }
        
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
          fprintf(file7, "%1.15e ", b);
          b = dataout[i] / (fabs(rfi) * rms);
          fprintf(file3, "%1.15e ", rms);
          fprintf(file5, "%1.15e ", av);
          

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
        fprintf(file3, "\n");
        fprintf(file5, "\n");
        fprintf(file7, "\n");
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
        printf("%d rms %f zwt %d / %d (or %d?) res %f %f z %f %f rms %f %f\n", j, rms, zwt, np, m, resmin, resmax, zmin, zmax, rmsmin, rmsmax);
        
        for(i=0; i<np; i++){
          fprintf(file2, "%f ", wtt[i]);
        }
        fprintf(file2, "\n");
      }
      fclose(file1);
      fclose(file2);
      fclose(file3);
      fclose(file5);
      fclose(file6);
      fclose(file7);
      
    }
    tav = polyfit(pfit, np, data, mcalc, wtt, dataout, fmode);
  } else {
    if (ns <= 0) {
      if (fabs(rfi) > 0.0) {
        for (j = 0; j < 10; j++) {
          rmsfilt(np, 16, fabs(rfi), data, wtt);  // was 8
          rms = rmscalc(np, data, wtt);
          //              for(i=0;i<np;i++) if(data[i] > fabs(rfi)*rms) wtt[i]=0;
          m = 0;
          for (i = 0; i < np; i++) {
            if (wtt[i] == 0) m++;
          }
          printf("%d rms %f zwt %d %d\n", j, rms, m, np);
        }
      }
    } else
      for (i = 0; i < np; i++) wtt[i] = wtts[i];
  }

  if (smooth) dsmooth(np, smooth, data, wtt);
  if (smooth)
    resol = (fstop * 1e06 / nspec) * smooth;
  else
    resol = (fstop * 1e06 / nspec);
  rms = rmscalc(np, data, wtt);
  //    rmstheory = 300.0*2.0/(sqrt(nline*(resol)*80.0*2.097e6*2.38e-9));  // 10 npci
  secint = nline * (double)nblock * (double)nspec * 2.0 / 420.0e6;
  // more complete expresssion for noise
  // include windowing
  printf("rmss %f theory %f np %d nline %d avpkpwr %f fracint %f fracfreq %f tav %f\n", rmss, rmstheory / (sqrt(resol * secint * 0.5)), np, nline,
         avpkpwr / nline, 3.0 * (line / 3.0 - nline) / line, (double)zwt / nnp, tav);
  if (np > 0) {
    //            if(smooth) dsmooth(np,smooth,dataout,wtt);
    if (pfit)
      plotfspec(np, fstart, dataout, wtt, water);
    else
      plotfspec(np, fstart, data, wtt, water);
  }
  if (water) {
    waterspec(2, np, data, 0, 0, 0, 0, wscalemax, wscalemin);
    waterwrite(2, np, data, yr, dy, hr, mn, sc, lst);
  }
  max = -1e99;
  min = 1e99;
  for (i = 0; i < np; i++)
    if (data[i] > max) max = data[i];
  for (i = 0; i < np; i++)
    if (data[i] < min) min = data[i];
  printf("max %f min %f maxp0 %6.3f minp0 %6.3f avp0 %6.3f\n", max, min, maxp0, minp0, avp0 / np0);
  toyrday(avsec / nline, &yr, &day, &hr, &mn, &sc);

  if (dgha != 48) {
    secs = avsec / nline;
    lst = (gst(secs) + lon) * 12.0 / PI;
    a = lst - (17.0 + 45.0 / 60.0 + 40.04 / 3600.0);
    secs = avsec / nline + (ghav / nline - a) * 3600.0;
    if (secs - avsec / nline > 43200) secs += -86400;
    if (secs - avsec / nline < -43200) secs += 86400;
    toyrday(secs, &yr, &day, &hr, &mn, &sc);
  }
  if (nline) specout(np, data, wtt, fname, nline, yr, day, hr, mn, sc, gha, dgha, ghav / nline, ghamax, ghamin);
  return 0;
}

double tempmon(double fstart, double fstep, double data[]) {
  int j, j1, j2, j3, j4, maxj;
  double freq, a, b, pwr, f1, f2, max, v, temp;
  j1 = j3 = (0.8 - fstart) / fstep + 0.5;
  j2 = j4 = (3.072 - fstart) / fstep + 0.5;
  max = -1e99;
  maxj = j1;
  for (j = j1; j < j2; j++) {
    j3 = j - 5;
    j4 = j + 5;
    if (j3 < j1) j3 = j1;
    if (j4 > j2) j4 = j2;
    a = data[j3];
    b = data[j4];
    f1 = fstart + j3 * fstep;
    f2 = fstart + j4 * fstep;
    freq = fstart + j * fstep;
    pwr = data[j] - a - (b - a) * (freq - f1) / (f2 - f1);
    if (pwr > max) {
      max = pwr;
      maxj = j;
    }
  }
  freq = fstart + maxj * fstep;
  v = (freq / 6.144) * 5.0;  // 2.5 volt ref
  temp = (v * 3.3 / 4.5 - 0.25) / 28e-3;
  temp = (freq / 0.884736 - 1.375) / 22.5e-3;
  temp = ((freq - 0.3072) / 0.7864 - 1.375) / 22.5e-3;
  //  printf("freq %f %f\n",freq,temp);
  return temp;
}

FILE *filew;
double secstart, xx2prev, hdone, yy1p;
void waterspec(int mode, int np, double data[], int dy, int hr, int mn, int sc, double wscalemax, double wscalemin)
// plot the spectrum
{
  char txt[256];

  int j, j1, j2, i, kk, jstep, njstep, njstep2;
  double dmax, dmin, dd, f, h, b;
  double xoffset, sum, sum2, sum3, hav;
  double wxscale, wyscale, woffset, ffstart, ffstop;
  double a, sec1, xx1, xx2, yy1, yy2, totp, ntot, ttotp, nntot;
  time_t now;

  if (mode == 0) {
    secstart = -1;
    hdone = 0;
    if ((filew = fopen("water.pos", "w")) == NULL) {
      printf("cannot open spe.pos:\n");
      return;
    }
    fprintf(filew, "%%!PS-Adobe-\n%c%cBoundingBox:  0 0 612 792\n%c%cEndProlog\n", '%', '%', '%', '%');
    fprintf(filew, "/cx1 {0 rlineto\n} def\n");
    fprintf(filew, "/cx2 {0 rlineto\n closepath\n } def\n");
    fprintf(filew, "/cy {1 1 sethsbcolor fill\n } def\n");
    fprintf(filew, "1 setlinewidth\n");
    //    fprintf(file, "2 setlinewidth\n");
    fprintf(filew, "/Times-Roman findfont\n 12 scalefont\n setfont\n");
  }
  xoffset = 80.0;
  dmax = -1.0e99;
  dmin = 1e99;
  dd = 0.0;
  j1 = 0;
  j2 = np;
  for (j = j1; j < j2; j++) {
    dd = data[j];
    if (dd > dmax) { dmax = dd; }
    if (dd < dmin) dmin = dd;
  }
  //    printf("dmax %f dmin %f\n",dmax,dmin);
  //   if(scaledb == 0.0 && dmax < 0.5) {dmax = 0.5; dmin = -0.5;}
  //   if(scaledb == 0.0 && dmax > 2000.0) dmax=2000.0;
  if (scaledb > 1.0) {
    i = (int)(dmin * 0.1);
    dmin = i * 10;
    dmin -= 10;
    dmax = dmin + scaledb;
  }
  //    printf("dmax %f dmin %f\n",dmax,dmin);
  if (mode == 1) {
    if (freqstart > 0.0 && freqstart < 80.0)
      ffstart = freqstart;
    else
      ffstart = 80.0;
    if (freqstop > 210.0)
      ffstop = 160.0;
    else
      ffstop = freqstop;
    ffstart = freqstart;
    ffstop = freqstop;
    wyscale = 0.25 * 2000.0 / (ffstop - ffstart);
    wxscale = 500.0;  // 512.0
    if (tstart >= 0 && tstop < 86400) wxscale = 500.0 * 24.0 / (tstop + 1 - tstart);
    woffset = 100.0;
    sec1 = dy * 24.0 * 3600.0 + hr * 3600.0 + mn * 60.0 + sc;
    if (secstart < 0.0) secstart = sec1;
    if (tstart >= 0 && tstop < 86400)
      xx2 = (sec1 - secstart) * wxscale / (24.0 * 3600.0) + 55.0;
    else
      xx2 = sec1 * wxscale / (24.0 * 3600.0) + 55.0;
    if (sec1 == secstart) {
      yy1p = woffset;
      xx1 = xx2 - 2;
      xx2prev = xx2;
      if ((ffstop - ffstart) > 200.0)
        b = 50.0;
      else
        b = 5.0;
      if ((fstop - ffstart) < 20.0) b = 1.0;
      for (a = 50; a <= 210; a += b) {
        if (a >= ffstart && a <= ffstop) {
          yy1 = (a - ffstart) * wyscale + woffset;
          fprintf(filew, "0 0 0 sethsbcolor\n");
          fprintf(filew, "newpath\n %6.3f %6.3f moveto \n %6.3f %6.3f lineto\n stroke\n", xx1 - 5.0, yy1, xx1, yy1);
          fprintf(filew, "%6.3f %6.3f moveto\n (%4.0f) show\n", xx1 - 31.0, yy1 - 4.0, a);
        }
      }

      a = (ffstop - ffstart) * wyscale + woffset + 100.0;
      yy1 = a;
      fprintf(filew, "0 0 0 sethsbcolor\n");
      fprintf(filew, "newpath\n %6.3f %6.3f moveto \n %6.3f %6.3f lineto\n stroke\n", xx1 - 5.0, yy1, xx1, yy1);
      fprintf(filew, "newpath\n %6.3f %6.3f moveto \n %6.3f %6.3f lineto\n stroke\n", xx1, yy1 - 100.0, xx1, yy1);
      fprintf(filew, "%6.3f %6.3f moveto\n (%4.0f) show\n", xx1 - 31.0, yy1 - 4.0, (a - ((ffstop - ffstart) * wyscale + woffset)) * 10.0);

      yy1 = 150.0 + woffset;
      fprintf(filew, "0 0 0 sethsbcolor\n");
      fprintf(filew, "%6.3f %6.3f moveto\n 90 rotate\n (Frequency (MHz)) show\n -90 rotate\n", xx1 - 33.0, yy1);

      for (j = 0; j <= 2000; j++) {
        yy1 = (j + 0) * 0.25 + woffset;
        yy2 = (j + 1) * 0.25 + woffset;
        h = (j / 2000.0) * 0.67;
        if (h > 0.67) h = 0.67;
        if (h < 0.0) h = 0.0;
        h = 0.67 - h;  // h=0 is red
        fprintf(filew, "newpath\n %6.2f %6.2f moveto\n", 568.0, yy1);
        fprintf(filew, "%6.3f 0 rlineto\n 0 %6.3f rlineto\n %6.3f 0 rlineto\n closepath\n", 5.0, yy2 - yy1, -5.0);
        fprintf(filew, "%5.3f 1 1 sethsbcolor fill\n", h);
      }
      for (j = 0; j <= 10; j++) {
        fprintf(filew, "0 0 0 sethsbcolor\n");
        yy1 = j * 50.0 + woffset;
        fprintf(filew, "newpath\n %6.3f %6.3f moveto \n %6.3f %6.3f lineto\n stroke\n", 568.0, yy1, 572.0, yy1);
        fprintf(filew, "/Times-Roman findfont\n 10 scalefont\n setfont\n");
        fprintf(filew, "%6.3f %6.3f moveto\n (%4.0f) show\n", 574.0, yy1 - 2.0, j * (wscalemax - wscalemin) / 10.0 + wscalemin);
        //                fprintf(filew, "%6.3f %6.3f moveto\n (10) show\n", 574.0, yy1 - 2.0);
        //                fprintf(filew, "/Times-Roman findfont\n 8 scalefont\n setfont\n");
        //                fprintf(filew, "%6.3f %6.3f moveto\n (%d) show\n", 583.0, yy1 + 4.0, j);
        fprintf(filew, "/Times-Roman findfont\n 12 scalefont\n setfont\n");
      }
      yy1 = 150.0 + woffset;
      //              fprintf(filew, "%6.3f %6.3f moveto\n 90 rotate\n (Delta temperture (K)) show\n -90 rotate\n", 595.0, yy1);
      fprintf(filew, "%6.3f %6.3f moveto\n 90 rotate\n (Temperature (K)) show\n -90 rotate\n", 595.0, yy1);
      fprintf(filew, "%6.3f %6.3f moveto\n 90 rotate\n (av. Ta (K)) show\n -90 rotate\n", xx1 - 33.0, 530.0 + woffset);
    }
    xx1 = xx2prev;
    xx2prev = xx2;
    if (smooth)
      jstep = smooth;
    else
      jstep = 1;
    totp = ntot = 0;
    ttotp = nntot = 0;
    for (j = 0; j < np; j += jstep) {
      f = j * fstep;
      if (f >= ffstart && f <= ffstop) {
        if (data[j] / ndataw < 2e03 * 0.5) {  // RFI effects total power
          totp += data[j] / ndataw;
          ntot++;
        }
        ttotp += data[j] / ndataw;
        nntot++;
        hav = sum3 = 0;
        for (kk = -jstep / 2; kk <= jstep / 2; kk++) {
          sum = sum2 = 0;
          njstep = 5.0 * 0.5 / fstep;  // average over 5 MHz
          njstep2 = njstep;
          if (njstep2 + j + kk >= np) njstep2 = np - j - kk;
          //                  for(k=-njstep;k<=njstep2;k++) { sum += data[j+kk+k]; sum2++; }
          //                  sum = sum/ndataw;
          //                  sum=sum/sum2;

          //                  h = (10.0*log10(data[j]/sum)) * 0.5;
          hav += data[j + kk] / ndataw - sum;
          sum3++;
        }
        h = hav / sum3;
        yy1 = (f - ffstart) * wyscale + woffset;
        yy2 = (f + fstep * jstep - ffstart) * wyscale + woffset;

        if (h > 0)
          h = ((hav / sum3) - wscalemin) * 0.67 / (wscalemax - wscalemin);  // wscale K full scale
        else
          h = 0;

        if (h > 0.67) h = 0.67;
        if (h < 0.0) h = 0.0;
        h = 0.67 - h;  // h=0 is red
        if (j == np - 10) printf("sec %f xx2 %f yy2 %f h %f num %d spec %f fstart %f fstop %f\n", sec1, xx2, yy2, h, np, f, ffstart, fstop);
        fprintf(filew, "newpath\n %6.2f %6.2f moveto\n", xx1, yy1);
        fprintf(filew, "%6.3f cx1\n 0 %6.3f rlineto\n %6.3f cx2\n", xx2 - xx1, yy2 - yy1, xx1 - xx2);
        fprintf(filew, "%5.3f cy\n", h);
      }
    }
    if (mn <= 30 && !hdone) {
      fprintf(filew, "0 0 0 sethsbcolor\n");
      fprintf(filew, "%6.3f %6.3f moveto\n (%02d) show\n", xx1 + (xx2 - xx1) * 0.5, woffset - 10.0, hr);
      hdone = 1;
    }
    if (mn > 30) hdone = 0;
    for (j = 0; j < np; j++) data[j] = 0;
    if (ntot > 0)
      totp = totp / ntot;
    else
      totp = ttotp / nntot;
    ndataw = 0;
    yy1 = (ffstop - ffstart) * wyscale + woffset + totp * 0.1;
    // printf("%02d:%02d:%02d totp %f\n",hr,mn,sc,totp);
    if (sec1 == secstart) yy1p = yy1;
    fprintf(filew, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xx1, yy1p, xx2, yy1);
    yy1p = yy1;
  }
  if (mode == 2) {
    fprintf(filew, "0 0 0 sethsbcolor\n");
    sprintf(txt, "file: %s", fname);
    fprintf(filew, "%6.2f %6.2f moveto\n (%s) show\n", 350.0 + xoffset, 50.0, txt);
    sprintf(txt, "UT %02d to %02d", tstart, tstop);
    fprintf(filew, "%6.2f %6.2f moveto\n (%s) show\n", 80.0 + xoffset, 50.0, txt);
    sprintf(txt, "fstart %3.0f fstop %3.0f pfit %d smooth %d resol %3.0f kHz rfi %3.1f nline %d secint %5.0f", freqstart, freqstop, pfit, smooth,
            resol * 1e-3, rfi, nline, secint);
    fprintf(filew, "%6.2f %6.2f moveto\n (%s) show\n", 0.0 + xoffset, 25.0, txt);
    now = time(NULL);
    fprintf(filew, "%d %d moveto\n (%s) show\n", 450, 35, ctime(&now));
    fprintf(filew, "showpage\n%c%cTrailer\n", '%', '%');
    fclose(filew);
  }
}

FILE *filewr;
void waterwrite(int mode, int np, double data[], int yr, int dy, int hr, int mn, int sc, double lst) {
  int i, j, jstep, k1, k2;
  double p, sum, num;

  if (mode == 0) {
    if ((filewr = fopen("water.txt", "w")) == NULL) {
      printf("cannot open water.txt:\n");
      return;
    }
  }
  if (mode == 1) {
    fprintf(filewr, "%4d:%03d:%02d:%02d:%02d       %7.4f", yr, dy, hr, mn, sc, lst);
    k1 = freqstart / (200.0 / np) + 0.5;
    k2 = freqstop / (200.0 / np) + 0.5;
    if (k2 > np) k2 = np;
    fprintf(filewr, " %7.0f %7.0f", freqstart, freqstop);
    if (smooth)
      jstep = smooth;
    else
      jstep = 1;
    i = 0;
    for (j = k1; j < k2; j += jstep) i++;
    fprintf(filewr, " %7d", i);
    for (j = k1; j < k2; j += jstep) {
      sum = num = 0;
      for (i = j - jstep / 2; i <= j + jstep / 2; i++) {
        if (i >= 0 && i < np) {
          sum += data[i];
          num++;
        }
      }
      if (num >= 0)
        p = sum / (ndataw * num);
      else
        p = 0;
      if (p <= 0) p = 1;
      fprintf(filewr, " %7.4f", 10.0 * log10(p));
    }
    fprintf(filewr, "\n");
  }
  if (mode == 2) { fclose(filewr); }
}

void plotfspec(int np, double fstart, double data[], double wtt[], int water)
// plot the spectrum
{
  char txt[256];
  int j, k, j1, j2, i, kk;
  double x, y, xp, yp, dmax, dmin, dd, f, fspan, totpp, scale, step, h, s, b;
  double xoffset, yoffset;
  time_t now;
  FILE *file;

  if ((file = fopen("spe.pos", "w")) == NULL) {
    printf("cannot open spe.pos:\n");
    return;
  }
  fprintf(file, "%%!PS-Adobe-\n%c%cBoundingBox:  0 0 612 700\n%c%cEndProlog\n", '%', '%', '%', '%');
  fprintf(file, "1 setlinewidth\n");
  //    fprintf(file, "2 setlinewidth\n");
  fprintf(file, "/Times-Roman findfont\n 12 scalefont\n setfont\n");

  xoffset = 80.0;
  yoffset = 100.0;
  fspan = np * fstep;
  dmax = -1.0e99;
  dmin = 1e99;
  dd = 0.0;
  j1 = 0;
  j2 = np;
  for (j = j1; j < j2; j++) {
    dd = data[j];
    if (dd > dmax && wtt[j] > 0.0) { dmax = dd; }
    if (dd < dmin && wtt[j] > 0.0) dmin = dd;
  }
  printf("dmax %f dmin %f\n", dmax, dmin);
  //   if(scaledb == 0.0 && dmax < 0.5) {dmax = 0.5; dmin = -0.5;}
  //   if(scaledb == 0.0 && dmax > 2000.0) dmax=2000.0;
  if (water) dmin = 0;
  if (scaledb > 1.0) {
    i = (int)(dmin * 0.1);
    dmin = i * 10;
    dmin -= 10;
    dmax = dmin + scaledb;
    scale = scaledb;
  } else {
    scale = dmax - dmin;
    dmin -= scale * 0.5;
    dmax += scale * 0.5;
    scale = dmax - dmin;
  }
  //    printf("dmax %f dmin %f\n",dmax,dmin);
  if (maxscale || minscale) {
    dmax = maxscale;
    dmin = minscale;
    scale = dmax - dmin;
  }
  for (y = 0; y < 2; y++)
    fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset, y * 480 + yoffset, xoffset + 400.0,
            y * 480 + yoffset);
  fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset, yoffset, xoffset, 480.0 + yoffset);
  fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset + 400.0, yoffset, xoffset + 400.0,
          480.0 + yoffset);
  fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xoffset + 0.0, yoffset, xoffset + 0.0,
          480.0 + yoffset);
  yp = 0;
  xp = 0;
  totpp = 0;
  kk = 0;
  for (k = 0; k < np; k++) {
    h = s = b = 0;
    x = k * 400.0 / (double)np;
    if (scale > 0.0) { totpp = data[k]; }
    if (wtt[k] == 0.0) {
      h = 0.6;
      s = 1;
      b = 1;
    }
    if (k > 0 && wtt[k - 1] == 0.0) {
      h = 0.6;
      s = 1;
      b = 1;
    }
    totpp = (totpp - dmin) / scale;
    y = totpp * 480.0;
    if (fstart + k * fstep > freqstart && fstart + k * fstep < freqstop)
      kk = 1;
    else
      kk = 0;
    if (y < 0) y = 0;
    if (y > 480) y = 480;
    if (k == 0) yp = y;
    if (y != yp) {
      if (kk)
        fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n %5.2f %5.2f %5.2f sethsbcolor stroke\n", xp + xoffset, yp + yoffset,
                x + xoffset, yp + yoffset, h, s, b);
      xp = x;
      if (y > yp) {
        if (kk)
          fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n  %5.2f %5.2f %5.2f sethsbcolor stroke\n", x + xoffset, yp + yoffset,
                  x + xoffset, y + yoffset, h, s, b);
      }
      if (yp > y) {
        if (kk)
          fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n %5.2f %5.2f %5.2f sethsbcolor stroke\n", x + xoffset, y + yoffset,
                  x + xoffset, yp + yoffset, h, s, b);
      }
    }
    yp = y;
  }
  if (kk)
    fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", xp + xoffset, yp + yoffset, 400.0 + xoffset,
            yp + yoffset);
  step = 20.0;
  if (fspan < 120.0) step = 10.0;
  if (fspan < 20.0) step = 2.0;
  if (fspan < 3.0) step = 0.2;
  if (fspan < 2.0) step = 0.1;
  //    i = (int) (fstart / fstep + 1.0);
  for (f = fstart; f <= fstart + fspan; f += step) {
    x = (f - fstart) * 400.0 / fspan;
    y = 0.0;
    fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, y + yoffset, x + xoffset,
            y - 10.0 + yoffset);
    if (step == 10.0)
      sprintf(txt, "%5.1f", f);
    else
      sprintf(txt, "%6.2f", f);
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 15.0, y - 20.0 + yoffset, txt);
  }
  if (scaledb > 0.0) {
    for (f = 0; f <= scaledb; f += 10) {
      x = 0;
      y = f * 480.0 / scaledb;
      fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, y + yoffset, x + xoffset + 5,
              y + yoffset);
      sprintf(txt, "%3.0f dB", dmin + f);
      fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 45.0, y + 1.0 + yoffset, txt);
    }
  } else {
    for (f = dmin; f < dmax; f += (dmax - dmin) * 0.1) {
      x = 0;
      y = (f - dmin) * 480.0 / (dmax - dmin);
      fprintf(file, "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 sethsbcolor stroke\n", x + xoffset, y + yoffset, x + xoffset + 5,
              y + yoffset);
      if (dmax - dmin < 10.0) {
        if (dmin > 300)
          sprintf(txt, "%4.2fK", f);
        else
          sprintf(txt, "%4.0f mK", f * 1e03);
      } else
        sprintf(txt, "%4.0f K", f);
      fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 50.0, y - 1.0 + yoffset, txt);
    }
  }
  sprintf(txt, "freq(MHZ)");
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 180.0, 65.0, txt);
  sprintf(txt, "file: %s", fname);
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", 350.0 + xoffset, 65.0, txt);
  sprintf(txt, "UT %02d to %02d", tstart, tstop);
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", 80.0 + xoffset, 65.0, txt);
  sprintf(txt, "fstart %3.0f fstop %3.0f pfit %d smooth %d resol %3.0f kHz rfi %3.1f int %d %5.0f sec rms %5.3f", freqstart, freqstop, pfit, smooth,
          resol * 1e-3, rfi, nline, secint, rmscalc(np, data, wtt));
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", 0.0 + xoffset, 50.0, txt);
  sprintf(txt, "peakpwr %5.3e pkpwrm %5.3e", peakpwr, pkpwrm);
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", 0.0 + xoffset, 35.0, txt);
  now = time(NULL);
  fprintf(file, "%d %d moveto\n (%s) show\n", 450, 20, ctime(&now));
  fprintf(file, "showpage\n%c%cTrailer\n", '%', '%');
  fclose(file);
}

double rmscalc(int np, double data[], double wt[]) {
  double a, b, c;
  int i;
  a = b = c = 0;
  for (i = 0; i < np; i++) {
    if (wt[i] > 0) {
      a += data[i] * wt[i];
      b += data[i] * data[i] * wt[i];
      c += wt[i];
    }
  }
  if (c > 0) {
    a = a / c;  // weighted mean
    if (b - a * a * c > 0.0)
      return sqrt((b - a * a * c) / c);  // assumes w = 0 or 1
    else {
      printf("rmserr %e a %e b %e c %e np %d\n", b - a * a * c, a, b, c, np);
      return 1e99;
    }
  } else
    return 0;
}
/*
void dsmooth(int np, int n, double data[], double wt[])
{
   double a,c,dd;
   double a00,a01,a10,a11,b0,b1,d,aa00,aa01,fre,w;
   int i,j,k;
   dd = 2.0 / n;
   for(i=0;i<np;i++){
   a00=a01=a10=a11=b0=b1=0;
   for(j=0;j<8*n+1;j++) {
     k = i + j - 4*n;
     if(k>=0 && k<np){
     fre = (k-i)*dd;
     c = exp(-fre*fre*0.69);
     w = c*wt[k];
     a00 += w;
     a11 += w * fre * fre;
     a01 += w * fre;
     b0 += w * data[k];
     b1 += w * data[k] * fre;
     }
    }
    a10 = a01;
    d = a00*a11 - a10*a01;
    aa00 =  a11/d;
//    aa10 = -a01/d;
    aa01 = -a10/d;
//    aa11 =  a00/d;
    a = aa00*b0 + aa01*b1;   // constant
//    b = aa10*b0 + aa11*b1;   // slope
    if(a00 > n/4.0) mcalc[i] = a;    // prevents smoothed freq with only a few freqs  n/4 best?
    else {mcalc[i] = data[i]; wt[i]=0;}
    }
   for(i=0;i<np;i++){
    if(wt[i] > 0.0)
       data[i]=mcalc[i];
   }
 }
*/

void dsmooth(int np, int n, double data[], double wt[]) {
  FILE *smoothed_weights = fopen("smoothed_weights.txt", "w");
  FILE *smoothed_weights_thresh = fopen("smoothed_weights_thresh.txt", "w");
  FILE *smoothed_weights_dec = fopen("smoothed_weights_decimated.txt", "w");

  double a, b, c, dd, fre, w;
  int i, j, k;
  dd = 2.0 / n;
  for (i = 0; i < np; i++) mcalc[i + np] = wt[i];
  for (i = 0; i < np; i++) {
    a = b = 0;
    for (j = 0; j < 8 * n + 1; j++) {
      k = i + j - 4 * n;
      if (k >= 0 && k < np) {
        fre = (k - i) * dd;
        c = exp(-fre * fre * 0.69);
        w = c * wt[k];
        a += w;
        b += w * data[k];
      }
    }
    fprintf(smoothed_weights, "%e\n",a);

    if (a > n / 4.0){
      mcalc[i] = b / a;  // prevents smoothed freq with only a few freqs  n/4 best?
      fprintf(smoothed_weights_thresh, "%e\n",a);
    } else {
      mcalc[i] = data[i];
      mcalc[i + np] = 0;
      fprintf(smoothed_weights_thresh, "%e\n",0.0);
    }
  }
  for (i = 0; i < np; i++) {
    if (wt[i] > 0.0) data[i] = mcalc[i];
    wt[i] = mcalc[i + np];
  }

  for (j = 0; j < np; j += smooth) {    
    fprintf(smoothed_weights_dec, "%e\n", wt[j]);
  }

}

double polyfit(int npoly, int nfreq, double ddata[], double mcalc[], double wtt[], double dataout[], int mode) {
  int i, j, k, kk, m1, m2;
  double re, max, dd;
  static long double aarr[100000], arrr[100000];
  static double bbrr[10000];
  short write = 1;

  for (i = 0; i < nfreq; i++) {
    kk = i * npoly;
    for (j = 0; j < npoly; j++) {
      mcalc[kk] = fitfun(i, j, mode);
      kk++;
    }
    if (wtt[i]==0) write = 0;
  }

  for (i = 0; i < npoly; i++) {
    re = 0.0;
    m1 = i;
    for (k = 0; k < nfreq; k++) {
      //            if(ddata[k] < 1) ddata[k]=1;
      //            dd = log10(ddata[k]);    // not clear how much log compression helps
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

      //     if(i==j) re += 1e-20;     // matrix inversion problem
      aarr[k] = arrr[k] = re;
    }
  }
  qrd(aarr, npoly, bbrr);
  //    MatrixInvert(npoly);
  max = 0;
  for (i = 0; i < npoly; i++) {
    for (j = 0; j < npoly; j++) {
      re = 0;
      for (m1 = 0; m1 < npoly; m1++) {
        k = m1 + i * npoly;
        kk = j + m1 * npoly;
        re += aarr[kk] * arrr[k];
      }
      //        if(i==j) printf("diag %d %e\n",i,re);
      //        if(i != j) printf("off %d %e\n",i,re);
      if (i == j)
        max += (re - 1.0) * (re - 1.0);
      else
        max += re * re;
    }
  }
  printf("maxtrix inv err %e\n", sqrt(max));
  for (i = 0; i < nfreq; i++) {
    re = 0.0;
    for (j = 0; j < npoly; j++) { re += bbrr[j] * fitfun(i, j, mode); }
    //        dd = pow(10.0,re);
    dd = re;
    dataout[i] = ddata[i] - dd;
    //        printf("ddata %f %f\n", ddata[i], re);
  }
  //    return pow(10.0,bbrr[0]);
  if(write){
    FILE *file = fopen("polyfit.txt", "w");
    for (i = 0; i < npoly; i++) {
      fprintf(file, "%g\n", bbrr[i]);
    }  
  }
  return bbrr[0];
}

double fitfun(int i, int j, int mode) {
  double a, f, f0, dfreq;
  f = fstart + i * fstep;
  f0 = fstart;
  dfreq = (freqstop - fstart) * 1.5;
  a = 1;
  if (j % 2) a = cos((f - f0) * 2.0 * PI * (j / 2 + 1) / dfreq);
  if (j % 2 == 0 && j > 0) a = sin((f - f0) * 2.0 * PI * (j / 2) / dfreq);
  if (mode == 1) a = pow(f / 100.0, -2.5 + j);
  if (mode == 2) a = pow(f, j);
  if (mode == 3) a = pow((f - (freqstop + fstart)) / 2.0, j);
  return a;
}

void qrd(long double a[], int n, double b[]) {
  int i, j, k;
  int pi, pk;
  long double c[NFIT], d[NFIT];
  long double scale, sigma, sum, tau;
  static long double qt[NFIT][NFIT], u[NFIT][NFIT];  // aa[30][30];

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

/*
void MatrixInvert(int nsiz)
{
    int ic0, id, i, j, ij, ic, n;
    long double re, mag, sumr;
    long double ttrr[1000];
    aarr[0] = 1.0 / aarr[0];
    bbrr[0] = bbrr[0] * aarr[0];

//     inversion by bordering
    for (n = 1; n <= nsiz - 1; n++) {
        ic0 = (n * (n + 1)) / 2;
        id = ic0 + n;
        for (i = 0; i < n; i++) {
            sumr = 0.0;
            ij = (i * (i + 1)) / 2; // A(I,0)
            for (j = 0; j < n; j++) {
                ic = ic0 + j;
                sumr += aarr[ic] * aarr[ij];
                if (j < i)
                    ij++;
                if (j >= i)
                    ij += j + 1;
            }
            ttrr[i] = sumr;
        }
        sumr = 0.0;
        for (i = 0; i < n; i++) {
            ic = ic0 + i;
            sumr += ttrr[i] * aarr[ic];
        }
        re = aarr[id] - sumr;
        mag = re * re;
        if (mag > 0.0)
            aarr[id] = re / mag;

        else
            aarr[id] = 0.0;
        sumr = 0.0;
        for (i = 0; i < n; i++) {
            ic = ic0 + i;
            sumr += aarr[ic] * bbrr[i];
        }
        bbrr[n] = aarr[id] * (bbrr[n] - sumr);
        for (i = 0; i < n; i++)
            bbrr[i] += -ttrr[i] * bbrr[n];
        for (i = 0; i < n; i++) {
            ic = ic0 + i;
            aarr[ic] = -aarr[id] * ttrr[i];
        }
        ij = 0;
        for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++) {
                ic = ic0 + j;
                aarr[ij] += -ttrr[i] * aarr[ic];
                ij++;
            }
        }
    }
}
*/

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

/* Convert Seconds to Yr/Day/Hr/Min/Sec */
void toyrday(double secs, int *pyear, int *pday, int *phr, int *pmin, int *psec) {
  double days, day, sec;
  int i;

  day = floor(secs / 86400.0);
  sec = secs - day * 86400.0;
  for (i = 1970; day > 365; i++) {
    days = ((i % 4 == 0 && i % 100 != 0) || i % 400 == 0) ? 366.0 : 365.0;
    day -= days;
  }
  *phr = sec / 3600.0;
  sec -= *phr * 3600.0;
  *pmin = sec / 60.0;
  *psec = sec - *pmin * 60;
  *pyear = i;
  day = day + 1;
  *pday = day;
  if (day == 366)  // fix for problem with day 366
  {
    days = ((i % 4 == 0 && i % 100 != 0) || i % 400 == 0) ? 366 : 365;
    if (days == 365) {
      day -= 365;
      *pday = day;
      *pyear = i + 1;
    }
  }
}

void rmsinit(double fun[]) {
  int i;
  for (i = 0; i < 32768; i++) fun[i] = pow((200.0 * i / 32768.0), -2.5);
}

double rmsf(double data[]) {
  int i, n1, n2;
  double sum, sum2, rms, x;
  sum = sum2 = rms = 0;
  // n1 = 60.0*32768.0/200.0; n2 = 80.0*32768.0/200.0;
  n1 = 60.0 * 32768.0 / 200.0;
  n2 = 80.0 * 32768.0 / 200.0;
  // n1 = 81.5*32768.0/200.0; n2 = 95.5*32768.0/200.0;
  for (i = n1; i < n2; i++) {
    sum += ffun[i] * data[i];
    sum2 += ffun[i] * ffun[i];
  }
  x = sum / sum2;
  for (i = n1; i < n2; i++) rms += (data[i] - x * ffun[i]) * (data[i] - x * ffun[i]);
  return sqrt(rms / (n2 - n1));
}

void rmsfilt(int np, int span, double level, double data[], double wt[]) {
  double a, b, c, rms;
  int i, j;
  for (j = 0; j < np - span; j++) {
    a = b = c = rms = 0;
    for (i = j; i < j + span; i++) {
      a += data[i] * wt[i];
      b += data[i] * data[i] * wt[i];
      c += wt[i];
    }
    if (c > 0) {
      a = a / c;                        // weighted mean
      rms = sqrt((b - a * a * c) / c);  // assumes w = 0 or 1
      for (i = j; i < j + span; i++) {
        if (data[i] - a > rms * level) wt[i] = 0;
      }
    }
  }
}

void b64init(int b64[]) {
  int i;
  for (i = 0; i < 256; i++) {
    b64[i] = 0;
    if (i >= 'A' && i <= 'Z') b64[i] = i - 'A';
    if (i >= 'a' && i <= 'z') b64[i] = i - 'a' + 26;
    if (i >= '0' && i <= '9') b64[i] = i - '0' + 52;
    if (i == '+') b64[i] = 62;
    if (i == '/') b64[i] = 63;
  }
}

void specout(int np, double data[], double wtt[], char fname[], int nline, int yr, int day, int hr, int min, int sc, double gha, double dgha,
             double ghav, double ghamax, double ghamin)
// plot the spectrum
{
  int j, k, m, jstep;
  double f;
  double fmpwr, a1, a2, av1, av2;
  FILE *file;

  if ((file = fopen("spe.txt", "w")) == NULL) {
    printf("cannot open spe.txt:\n");
    return;
  }
  av1 = av2 = 0;
  a1 = a2 = 1e-99;
  for (j = 0; j < np; j++) {
    f = fstart + j * fstep;
    if (f > 87.5 && f < 108 && wtt[j] == 0) {
      av1 += data[j];
      a1++;
      m = 1;
      for (k = j + 1; k < np && m; k++) {
        if (wtt[k]) {
          av2 += data[k];
          a2++;
          m = 0;
        }
      }
      m = 1;
      for (k = j - 1; k >= 0 && m; k--) {
        if (wtt[k]) {
          av2 += data[k];
          a2++;
          m = 0;
        }
      }
    }
  }
  fmpwr = av1 / a1 - av2 / a2;
  if (smooth)
    jstep = smooth;
  else
    jstep = 1;
  for (j = 0; j < np; j += jstep) {
    f = fstart + j * fstep;
    if (j == 0)
      fprintf(file, "%12.6f %12.6f %4.0f %d // %s %04d:%03d:%02d:%02d:%02d %7.3f %7.3f %7.3f %7.3f %7.3f fmpwr %7.3f\n", f, data[j], wtt[j], nline,
              fname, yr, day, hr, min, sc, gha, dgha, ghav, ghamax, ghamin, fmpwr);
    else
      fprintf(file, "%12.6f %12.6f %4.0f\n", f, data[j], wtt[j]);
  }
  fclose(file);
}

double gst(double ttime) {
  double secs, pdum;
  int i;
  secs = (1999 - 1970) * 31536000.0 + 17.0 * 3600.0 + 16.0 * 60.0 + 20.1948;
  for (i = 1970; i < 1999; i++) {
    if ((i % 4 == 0 && i % 100 != 0) || i % 400 == 0) secs += 86400.0;
  }

  return (modf((ttime - secs) / 86164.09053, &pdum) * TWOPI);
  /* 17 16 20.1948 UT at 0hr newyear1999 */
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

void radec_azel(double ha, double dec, double latt, double *azs, double *elevs) {
  /* convert from sky to antenna coords (azel mount) */
  /* input: ha,dec,latt
     output: azs=azimuth of source
     elevs=elevation of source
   */
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
