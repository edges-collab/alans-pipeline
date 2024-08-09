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
#define NDATA 100000
#define NBEAM 256
#define NFIT 50  // more than 50 may not be accurate

void plotfspec(int, double *, double *, double *, double *, double, int, int, char *, char *, char *);
void outfspec(int, double *, double *, double *, double *, double, int, int, char *);
void outcal(int, double *, complex double *, double *, double *, double *, double *, double *, double *, char *);
void outsim(int, double *, double *, int, int, int, int, int);
void readcal(int *, double *, complex double *, double *, double *, double *, double *, double *, double *);
void plotnwave(int, double *, double *, double *, double *, double, int);
double rmscalc(int, double *, double *, double *);
void dsmooth(int, double *, double *, double *, double *, int);
double polyfitr(int, int, double *, double *, double *, double *, double *, long double *, double *);
double polyfitrc(int, int, double *, double *, double *, double *, double *, long double *, double *, int);
double polyfitr2(int, int, double *, double *, double *, double *, double *, long double *, double *, double *);
double fitfun(int, int, int, double *);
double wavemodel(complex double, complex double, double, double, double, double, double);
double w3p(complex double, complex double, double, double, double, double, double, double);
double w3pinv(complex double, complex double, double, double, double, double, double, double);
void fittp(int, double *, complex double *, double *, int, double *, complex double *, int, double *, double *, long double *, double *, int, int,
           char *);
void plotvna(int, double *, complex double *, double *, int, double *, complex double *, double *, double *, int, double, double, int, char *);
void wavefit(int, int, int, double *, double *, double *, complex double *, complex double *, double, double, double, double *, double *, double *,
             double *, double *, double *, long double *, double *, double *);
double loss(complex double, double);
double lossmodel(double, complex double, int);
double lossmodel2(double, complex double);
double rigloss(complex double, complex double, complex double, complex double);
complex double cabl2(double, double, complex double, int, complex double *, complex double *, double);
double antloss(double, double *, double *, long double *, double *, double, int);
double lossinv(double, double, double);
double beamcorr(double, int, double, double, double *, double *, long double *, double *, int, int, double *, double *, double *, int, double *,
                double, double, int, int);
double beamcorr2(double, int, double, double, double *, double *, long double *, double *, int, int, double *, double *, double *, int, double *,
                 double, double, int);
double gst(double);
void radec_azel(double, double, double, double *, double *);
void radec_azel2(double, double, double, double, double, double *, double *);
void GalactictoRadec(double, double, double *, double *);
void sunradec(double, double *, double *);
double tosecs(int, int, int, int, int);
void toyrday(double, int *, int *, int *, int *, int *);
void qrd(long double *, int, double *);
void pix2ang_ring(const long, long, double *, double *);
double gauss(void);

double parm1, parm2, parm3, parm4;

int main(int argc, char *argv[]) {
  static long double aarr[NDATA];
  static double bbrr[NDATA], mcalc[NDATA * NFIT];
  static double freqant[NDATA], freqopen[NDATA], freqshort[NDATA], freqhot[NDATA], freqamb[NDATA], freqs11ant[NDATA];
  static double freqs11lna[NDATA], freqs11hot[NDATA], freqs11amb[NDATA], freqs11cab1[NDATA], freqs11cab2[NDATA], freqs11rig[NDATA];
  static double wtant[NDATA], wtopen[NDATA], wtshort[NDATA], wthot[NDATA], wtamb[NDATA], wttant[NDATA], wttlna[NDATA], wtthot[NDATA], wttamb[NDATA],
      wttcab1[NDATA], wttcab2[NDATA], wttrig[NDATA];
  static double spant[NDATA], spopen[NDATA], spshort[NDATA], sphot[NDATA], spamb[NDATA], data[NDATA], dataout[NDATA], dataout2[NDATA],
      fitf[NDATA * NFIT], md[NDATA], ltemp[NDATA];
  static double sspant[NDATA], sspopen[NDATA], ssphot[NDATA], sspamb[NDATA],
      wtemp[NDATA];  // temporary - to preserve original
  static double tlna0[NDATA], tlna1[NDATA], tlna2[NDATA], skymodel[NDATA], freqcal[NDATA];
  static complex double s11ant[NDATA], s11lna[NDATA], s11hot[NDATA], s11amb[NDATA], s11cab1[NDATA], s11cab2[NDATA], s11rig[NDATA], s12rig[NDATA],
      s22rig[NDATA], T;
  static complex double ss11ant[NDATA], ss11lna[NDATA], ss11hot[NDATA], ss11amb[NDATA], ss11cab[NDATA], ss11rig[NDATA], ss12rig[NDATA],
      ss22rig[NDATA], tmp[NDATA], tmpsm[NDATA];
  static double bb[NBEAM], bb2[NBEAM];

  double freq, wt, re, re1, re2, re3, re4, re5, re6, re7, re8, re9, re10, re11, re12, re13, re14, im, im1, im2, im3, im4, im5, im6, im7, im8, im9,
      im10, im11, im12, im13, im14, bspac;
  double tt, La, Lh, Lhh, rms, rms1, rms2, rms3, fstart, fstop, wfstart, wfstop, fcen, eorcen, eorwid, eoramp, atten, dlyrob, days, adb, ldb, gha,
      dgha, specin, freqr, antaz, dscale;
  double spe_ind, spe_inderr, gamma, gammaerr, ion_abs, ion_em, ssun, sunind, ghav, ghamax, ghamin, ion, tau, noise, opn, fbstart, fbstop, lr, li, ar,
      ai;
  double thot, tamb, tcold, tant, tload, tcab, s, sca, s1, s2, ofs, err1, lst, secs, aloss, delaylna, delayln, delaycorr, dbcorr, t150, delayant,
      tant2, freqref;
  int i, j, k, kk, nant, nopen, nshort, nhot, namb, mode, ns11ant, ns11lna, ns11hot, ns11amb, ns11cab1, ns11cab2, ns11rig, nfit, nfit1, imp, db,
      nfit2, nfit3, nfit4, iter, mfit, cfit, wfit, skymode, yr, dy, hr, mn, sc, smooth, cal, wtmode, lmode, cmb, nline, site, map;
  int ncal, cons, binteg, nter, nocal, bfit, bfit2, mdd, skymod2, test, i1, i2, i3, i4, i5, i6, i7, yrh, dyh, hrh, mnh, sch, yrc, dyc, hrc, mnc, scc,
      yro, dyo, hro, mno, sco, yrs, dys, hrs, mns, scs, low, sim, rr, lna_poly;
  static double tcal_sca[NDATA], tcal_ofs[NDATA], tcal_scasm[NDATA], tcal_ofssm[NDATA], wtcal[NDATA];
  FILE *file1;
  char fname[256], buf[2048], title[256], hottitle[256], ambtitle[256], opentitle[256], shorttitle[256], datfile[256], ambfile[256], hotfile[256],
      openfile[256], shortfile[256], info[256];
  char antfname[256], hotfname[256], ambfname[256], openfname[256], shortfname[256], rigfname[256], lnafname[256];
  mode = 0;
  nant = nopen = nshort = nhot = namb = ns11ant = ns11lna = ns11hot = ns11amb = ns11cab1 = ns11cab2 = ns11rig = 0;
  wtmode = 0;
  cons = 0;
  nocal = 0;
  low = 0;
  sim = ion = tau = noise = rr = 0;
  lr = li = ar = ai = 0;
  datfile[0] = ambfile[0] = hotfile[0] = openfile[0] = shortfile[0] = 0;
  title[0] = hottitle[0] = ambtitle[0] = opentitle[0] = shorttitle[0] = 0;
  antaz = 90.0;
  tamb = 300.0;   // should agree with what is in acqplot2 - but doesn't really
                  // matter
  tload = 300.0;  // internal comparison load - value doesn't matter as long as
                  // it have same value as used in calibration
  cal = 0;
  tcab = 0;
  atten = 0;
  lna_poly = -1;
  Lhh = 1;
  wfstart = 50;
  wfstop = 200;
  imp = 0;
  db = 0;
  secs = 0;
  yr = yrh = yrc = yro = yrs = dy = hr = mn = sc = cmb = nline = 0;
  fstart = 50;
  fstop = 200;
  fbstart = 0;
  site = 0;
  map = 0;
  lst = 0;
  skymode = skymod2 = -1;
  tcold = tant = 273.0 + 25.0;
  aloss = 0;
  bfit = bfit2 = 0;
  mfit = 1;
  thot = 273.0 + 95.0;
  delaylna = 0;
  wfit = 7;
  cfit = 7;
  nfit1 = 27;
  nfit2 = 27;
  nfit3 = 27;
  nfit4 = 37;
  dlyrob = 0.5 / 3e8;
  smooth = 1;
  Lh = 0.05;  // estimated loss of 6 inches 0.141 stainless steel dB at 100 MHz
              // thot  and tcold have to be correct or both scaled 1K error ~
              // 10mK
  eorcen = eorwid = eoramp = 0;
  delaycorr = dbcorr = 0;
  lmode = -1;
  adb = ldb = opn = 0;
  binteg = 0;
  dscale = 0;
  nter = 8;
  delayant = 0;
  mdd = 0;
  test = 0;
  specin = -2.5;
  sunind = 0.5;
  for (i = 0; i < argc - 1; i++) {
    sscanf(argv[i], "%79s", buf);
    if (strstr(buf, "-fstart")) { sscanf(argv[i + 1], "%lf", &fstart); }  // need to be first is used
    if (strstr(buf, "-fstop")) { sscanf(argv[i + 1], "%lf", &fstop); }    // need to be first if used
    if (strstr(buf, "-spant")) mode = 1;
    if (strstr(buf, "-spopen")) mode = 2;
    if (strstr(buf, "-spshort")) mode = 5;
    if (strstr(buf, "-sphot")) mode = 3;
    if (strstr(buf, "-spcold")) mode = 4;
    if (strstr(buf, "-s11ant")) mode = 11;
    if (strstr(buf, "-s11open")) mode = 12;
    if (strstr(buf, "-s11short")) mode = 16;
    if (strstr(buf, "-s11lna")) mode = 13;
    if (strstr(buf, "-s11hot")) mode = 14;
    if (strstr(buf, "-s11cold")) mode = 15;
    if (strstr(buf, "-cals11")) mode = 17;
    if (strstr(buf, "-cals11_cold")) mode = 18;
    if (strstr(buf, "-cals11_hot")) mode = 19;
    if (strstr(buf, "-cals11_open")) mode = 20;
    if (strstr(buf, "-cals11_short")) mode = 21;
    if (strstr(buf, "-cals11_sim1")) mode = 22;
    if (strstr(buf, "-cals11_sim2")) mode = 23;
    if (strstr(buf, "-cals11_sim3")) mode = 24;
    if (strstr(buf, "-cals11_sim4")) mode = 25;
    if (strstr(buf, "-cals11_cbox")) mode = 26;
    if (strstr(buf, "-cals11_3db")) mode = 27;
    if (strstr(buf, "-cals11_cbx")) mode = 28;
    if (strstr(buf, "-cals11_noise")) mode = 29;
    if (strstr(buf, "-cals11_s11rig")) mode = 30;  // just get s11rig
    if (strstr(buf, "-hotsparms"))
      mode = 31;                                 // from /data5/edges/data/CalHotLoadCableData/ but has to be
                                                 // fit if freq spacing different
    if (strstr(buf, "-cals11_sm33")) mode = 33;  // more recent calibrations just using simulator 3
    if (strstr(buf, "-nocal")) nocal = 1;
    if (strstr(buf, "-ydhms")) { sscanf(argv[i + 1], "%d:%d:%d:%d:%d", &yr, &dy, &hr, &mn, &sc); }
    if (strstr(buf, "-skymode")) {
      sscanf(argv[i + 1], "%d", &skymode);
    }  // -1 = no beamcorrection 0 = beamcorrection 256 = use beamcorr as data +
       // 1,2,4,8,16,32,64,128 opts
    if (strstr(buf, "-skymod2")) sscanf(argv[i + 1], "%d",
                                        &skymod2);  // use for new slower beamcorrection using new Haslam map
    if (strstr(buf, "-lmode"))
      sscanf(argv[i + 1], "%d", &lmode);  // -1 = to use aloss 0=for hot load
                                          // 2=hot+rigloss 1= balun + bean using
                                          // interpolated antloss 5or6=lowband balun +conn
    if (strstr(buf, "-wtmode")) {
      sscanf(argv[i + 1], "%d", &wtmode);
    }                                                                        // sets weighting for fit and > 10 for  calibration 1 = normal zero wt
                                                                             // outside wfstart & wfstop
    if (strstr(buf, "-nfit1")) { sscanf(argv[i + 1], "%d", &nfit1); }        // fitting to calibration data in specal.txt  default=27
    if (strstr(buf, "-nfit2")) { sscanf(argv[i + 1], "%d", &nfit2); }        // calibration fitting of s11 for hot,cold and cable  default=27
    if (strstr(buf, "-nfit3")) { sscanf(argv[i + 1], "%d", &nfit3); }        // LNA s11 fitting default=27
    if (strstr(buf, "-lna_poly")) { sscanf(argv[i + 1], "%d", &lna_poly); }  // LNA s11 fitting default=27
    if (strstr(buf, "-nfit4")) { sscanf(argv[i + 1], "%d", &nfit4); }        // antenna s11 fitting default=37
    if (strstr(buf, "-mfit")) { sscanf(argv[i + 1], "%d", &mfit); }          // num parameters in fit to spectrum
    if (strstr(buf, "-wfit")) { sscanf(argv[i + 1], "%d", &wfit); }          // noisewave fitting default=4
    if (strstr(buf, "-bfit")) { sscanf(argv[i + 1], "%d", &bfit); }          // optional fitting of beam data default = 0
    if (strstr(buf, "-cfit")) { sscanf(argv[i + 1], "%d", &cfit); }          // final calibration fitting not used if specal.txt is used default=7
    if (strstr(buf, "-smooth")) { sscanf(argv[i + 1], "%d", &smooth); }      // optional smoothing of antenna spectrum after fitting
    if (strstr(buf, "-cons")) { sscanf(argv[i + 1], "%d", &cons); }          // optional fitting constaints
    if (strstr(buf, "-aloss")) { sscanf(argv[i + 1], "%lf", &aloss); }       // alternate simple antenna loss needs lmode = -1 to activate
    if (strstr(buf, "-atten")) { sscanf(argv[i + 1], "%lf", &atten); }       // attenuation between antenna and LNA
    if (strstr(buf, "-adb")) { sscanf(argv[i + 1], "%lf", &adb); }           // correction to ant s11 in dB
    if (strstr(buf, "-ldb")) { sscanf(argv[i + 1], "%lf", &ldb); }           // correction to lna s11 in dB
    if (strstr(buf, "-lr")) { sscanf(argv[i + 1], "%lf", &lr); }             // correction to lna s11 in re,im
    if (strstr(buf, "-li")) { sscanf(argv[i + 1], "%lf", &li); }             // correction to lna s11 in re,im
    if (strstr(buf, "-ar")) { sscanf(argv[i + 1], "%lf", &ar); }             // correction to ant s11 in re,im
    if (strstr(buf, "-ai")) { sscanf(argv[i + 1], "%lf", &ai); }             // correction to ant s11 in re,im
    if (strstr(buf, "-opn")) { sscanf(argv[i + 1], "%lf", &opn); }           // correction to lna to open at VNA
    if (strstr(buf, "-thot")) { sscanf(argv[i + 1], "%lf", &thot); }         // temperature of hot calibration load
    if (strstr(buf, "-tcold")) { sscanf(argv[i + 1], "%lf", &tcold); }       // temperature of cold (normally room temperature) calibration load
    if (strstr(buf, "-tcab")) {
      sscanf(argv[i + 1], "%lf", &tcab);
    }  // temperature of cable used for noise wave calibration get set to tcold
       // unless entered
    if (strstr(buf, "-tant")) {
      sscanf(argv[i + 1], "%lf", &tant);
    }                                                                           // temperature of antenna used in loss calculation needed even when
                                                                                // specal.txt is used
    if (strstr(buf, "-wfstart")) { sscanf(argv[i + 1], "%lf", &wfstart); }      // start of weighted data default=50
    if (strstr(buf, "-wfstop")) { sscanf(argv[i + 1], "%lf", &wfstop); }        // stop of weighted data default=200
    if (strstr(buf, "-fbstart")) { sscanf(argv[i + 1], "%lf", &fbstart); }      // beamfit start
    if (strstr(buf, "-fbstop")) { sscanf(argv[i + 1], "%lf", &fbstop); }        // beamfit stop
    if (strstr(buf, "-delaylna")) { sscanf(argv[i + 1], "%lf", &delaylna); }    // adapter correction on LNA S11
    if (strstr(buf, "-delayant")) { sscanf(argv[i + 1], "%lf", &delayant); }    // correction for antenna
    if (strstr(buf, "-dlyrob")) { sscanf(argv[i + 1], "%lf", &dlyrob); }        // one-way delay in Roberts balun
    if (strstr(buf, "-delaycorr")) { sscanf(argv[i + 1], "%lf", &delaycorr); }  // VNA corr
    if (strstr(buf, "-dbcorr")) { sscanf(argv[i + 1], "%lf", &dbcorr); }        // VNA corr
    if (strstr(buf, "-Lh")) {
      sscanf(argv[i + 1], "%lf", &Lh);
    }                                                                     // loss of hot load at 100 MHz uses tamb for loss calc Lh = -1 for hot
                                                                          // load model Lh = -2 to use s12,s22
    if (strstr(buf, "-eorcen")) { sscanf(argv[i + 1], "%lf", &eorcen); }  // eor model center freq 0=no EOR search forces different set of parms
    if (strstr(buf, "-eorwid")) { sscanf(argv[i + 1], "%lf", &eorwid); }  // eor model FWHM width 0=no EOR search
    if (strstr(buf, "-eoramp")) { sscanf(argv[i + 1], "%lf", &eoramp); }  // eor model FWHM width 0=no EOR search
    if (strstr(buf, "-antaz")) { sscanf(argv[i + 1], "%lf", &antaz); }    // antenna azimuth
    if (strstr(buf, "-dscale")) { sscanf(argv[i + 1], "%lf", &dscale); }  // fix plot scale
    if (strstr(buf, "-specin")) { sscanf(argv[i + 1], "%lf", &specin); }  // assumed foreground spectral index
    if (strstr(buf, "-sunind")) { sscanf(argv[i + 1], "%lf", &sunind); }  // assumed sun spectral index
    if (strstr(buf, "-binteg")) { sscanf(argv[i + 1], "%d", &binteg); }   // beamcorrection time span is secs default 1
    if (strstr(buf, "-mdd")) { sscanf(argv[i + 1], "%d", &mdd); }         // defalt=mdd=0 pow(f/150,-2.5+i)  mdd=1 for gamma etc.
    if (strstr(buf, "-sim")) { sscanf(argv[i + 1], "%d", &sim); }         // simulate data
    if (strstr(buf, "-test")) { sscanf(argv[i + 1], "%d", &test); }       // test = 1 use beamcorr as data
    if (strstr(buf, "-cmb")) { sscanf(argv[i + 1], "%d", &cmb); }         // subtraction of cmb for spectral index and beam correction
    if (strstr(buf, "-parm1")) { sscanf(argv[i + 1], "%lf", &parm1); }    // search for optimum values
    if (strstr(buf, "-parm2")) { sscanf(argv[i + 1], "%lf", &parm2); }    // search for optimum values
    if (strstr(buf, "-parm3")) { sscanf(argv[i + 1], "%lf", &parm3); }    // search for optimum values
    if (strstr(buf, "-parm4")) { sscanf(argv[i + 1], "%lf", &parm4); }    // search for optimum values
    if (strstr(buf, "-ion")) { sscanf(argv[i + 1], "%lf", &ion); }        // add ionosphere
    if (strstr(buf, "-tau")) { sscanf(argv[i + 1], "%lf", &tau); }        // opacity
    if (strstr(buf, "-noise")) { sscanf(argv[i + 1], "%lf", &noise); }    // add noise
    if (strstr(buf, "-rr")) {
      sscanf(argv[i + 1], "%d", &rr);
      if (rr) srand(rr);
    }                                                                // srand
    if (strstr(buf, "-site")) { sscanf(argv[i + 1], "%d", &site); }  // 0 = MRO 1 = Oregon
    if (strstr(buf, "-map")) { sscanf(argv[i + 1], "%d", &map); }    // 0 = Haslam
    if (mode) {
      sscanf(argv[i + 1], "%s", fname);
      printf("Reading %s\n", fname);

      if ((file1 = fopen(fname, "r")) == NULL) {
        printf("cannot open file:%s\n", fname);
        return 0;
      }
      k = j = imp = db = 0;
      while (fgets(buf, 2048, file1) != 0) {
        if (mode >= 1 && mode < 10) {
          sscanf(buf, "%lf %lf %lf", &freq, &tt, &wt);
          if (mode == 1 && j == 0) {
            sscanf(buf, "%*s %*s %*s %d %*s %s %d:%d:%d:%d:%d %lf %lf %lf %lf %lf", &nline, datfile, &yr, &dy, &hr, &mn, &sc, &gha, &dgha, &ghav,
                   &ghamax, &ghamin);
          }
          if (mode == 2 && j == 0) {
            sscanf(buf, "%*s %*s %*s %*s %*s %s %d:%d:%d:%d:%d %lf %lf %lf %lf %lf", openfile, &yro, &dyo, &hro, &mno, &sco, &gha, &dgha, &ghav,
                   &ghamax, &ghamin);
          }
          if (mode == 5 && j == 0) {
            sscanf(buf, "%*s %*s %*s %*s %*s %s %d:%d:%d:%d:%d %lf %lf %lf %lf %lf", shortfile, &yrs, &dys, &hrs, &mns, &scs, &gha, &dgha, &ghav,
                   &ghamax, &ghamin);
          }
          if (mode == 3 && j == 0) {
            sscanf(buf, "%*s %*s %*s %*s %*s %s %d:%d:%d:%d:%d %lf %lf %lf %lf %lf", hotfile, &yrh, &dyh, &hrh, &mnh, &sch, &gha, &dgha, &ghav,
                   &ghamax, &ghamin);
          }
          if (mode == 4 && j == 0) {
            sscanf(buf, "%*s %*s %*s %*s %*s %s %d:%d:%d:%d:%d %lf %lf %lf %lf %lf", ambfile, &yrc, &dyc, &hrc, &mnc, &scc, &gha, &dgha, &ghav,
                   &ghamax, &ghamin);
          }
          if (freq >= fstart && freq <= fstop) {
            if (mode == 1) {
              freqant[j] = freq;
              md[j] = spant[j] = tt;
              wtant[j] = wt;
              nant = j + 1;
            }
            if (mode == 2) {
              freqopen[j] = freq;
              spopen[j] = tt;
              wtopen[j] = wt;
              nopen = j + 1;
            }
            if (mode == 5) {
              freqshort[j] = freq;
              spshort[j] = tt;
              wtshort[j] = wt;
              nshort = j + 1;
            }
            if (mode == 3) {
              freqhot[j] = freq;
              sphot[j] = tt;
              wthot[j] = wt;
              nhot = j + 1;
            }
            if (mode == 4) {
              freqamb[j] = freq;
              spamb[j] = tt;
              wtamb[j] = wt;
              namb = j + 1;
            }
            j++;
          }
        }
        if (mode > 10 && mode < 17) {
          if (strstr(buf, "IMP")) imp = 1;
          if (strstr(buf, "DB")) db = 1;
          if (strstr(buf, "BEGIN")) k = 1;
          if (strstr(buf, "END")) k = 0;
          sscanf(buf, "%lf,%lf,%lf", &freq, &re, &im);
          if (k >= 2) {
            T = re + im * I;
            if (imp) T = (T - 50.0) / (T + 50.0);
            if (db) T = pow(10.0, 0.05 * re) * cexp(im * PI * I / 180.0);
            if (freq >= fstart * 1e6 && freq <= fstop * 1e6) {
              if (mode == 11) {
                s11ant[j] = T;
                freqs11ant[j] = freq / 1e6;
                wttant[j] = 1;
                ns11ant = j + 1;
                strcpy(antfname, fname);
              }
              if (mode == 12) {
                s11cab1[j] = T;
                freqs11cab1[j] = freq / 1e6;
                wttcab1[j] = 1;
                ns11cab1 = j + 1;
                strcpy(openfname, fname);
              }
              if (mode == 16) {
                s11cab2[j] = T;
                freqs11cab2[j] = freq / 1e6;
                wttcab2[j] = 1;
                ns11cab2 = j + 1;
                strcpy(shortfname, fname);
              }
              if (mode == 13) {
                s11lna[j] = T;
                freqs11lna[j] = freq / 1e6;
                wttlna[j] = 1;
                ns11lna = j + 1;
                strcpy(lnafname, fname);
              }
              if (mode == 14) {
                s11hot[j] = T;
                freqs11hot[j] = freq / 1e6;
                wtthot[j] = 1;
                ns11hot = j + 1;
                strcpy(hotfname, fname);
              }
              if (mode == 15) {
                s11amb[j] = T;
                freqs11amb[j] = freq / 1e6;
                wttamb[j] = 1;
                ns11amb = j + 1;
                strcpy(ambfname, fname);
              }
              if (buf[0] != 'B') j++;
            }
          }
          if (k) k++;
        }
        if (mode >= 17 && mode < 40) {  // Raul's format for all s11 cal file
          // freq [MHz],real(LNA),imag(LNA),real(ambient),imag(ambient),real(hot),
          // imag(hot),real(open_cable),imag(open_cable),real(shorted_cable),
          // imag(shorted_cable),real(s11_semirigid),imag(s11_semirigid),
          // real(s12s21_semirigid),imag(s12s21_semirigid),real(s22_semirigid),
          // imag(s22_semirigid),real(ant_sim1),imag(ant_sim1),real(ant_sim2),imag(ant_sim2)
          //
          // noise_source_s-parameters.txt
          // # freq [MHz]
          // re(noise_source_plus_3dB),im,re(noise_source_plus_3dB_plus_cable_box),im,re(s11_cable_box),im,re(s12s21_cable_box),im,re(s22_cable_box),im
          //
          freq = 0;  // ensure freq is read
          if (buf[0] != '#') {
            sscanf(buf,
                   "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
                   "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                   &freq, &re, &im, &re1, &im1, &re2, &im2, &re3, &im3, &re4, &im4, &re5, &im5, &re6, &im6, &re7, &im7, &re8, &im8, &re9, &im9, &re10,
                   &im10, &re11, &im11, &re12, &im12, &re13, &im13, &re14, &im14);
            if (freq >= fstart && freq <= fstop) {
              if (mode == 17) {  // don't overwrite other cals
                T = re + im * I;
                s11lna[j] = T;
                freqs11lna[j] = freq;
                wttlna[j] = 1;
                ns11lna = j + 1;
                strcpy(lnafname, fname);
                T = re1 + im1 * I;
                s11amb[j] = T;
                freqs11amb[j] = freq;
                wttamb[j] = 1;
                ns11amb = j + 1;
                strcpy(ambfname, fname);
                T = re2 + im2 * I;
                s11hot[j] = T;
                freqs11hot[j] = freq;
                wtthot[j] = 1;
                ns11hot = j + 1;
                strcpy(hotfname, fname);
                T = re3 + im3 * I;
                s11cab1[j] = T;
                freqs11cab1[j] = freq;
                wttcab1[j] = 1;
                ns11cab1 = j + 1;
                strcpy(openfname, fname);
                T = re4 + im4 * I;
                s11cab2[j] = T;
                freqs11cab2[j] = freq;
                wttcab2[j] = 1;
                ns11cab2 = j + 1;
                strcpy(shortfname, fname);
              }
              if (mode == 17 || mode == 30) {  // all calfile
                T = re5 + im5 * I;
                s11rig[j] = T;
                freqs11rig[j] = freq;
                wttrig[j] = 1;
                ns11rig = j + 1;
                strcpy(rigfname, fname);
                T = re6 + im6 * I;
                s12rig[j] = T;
                freqs11rig[j] = freq;
                wttrig[j] = 1;
                ns11rig = j + 1;
                strcpy(rigfname, fname);
                T = re7 + im7 * I;
                s22rig[j] = T;
                freqs11rig[j] = freq;
                wttrig[j] = 1;
                ns11rig = j + 1;
                strcpy(rigfname, fname);
              }
              if (mode > 17 && mode < 26) {
                if (mode == 18) T = re1 + im1 * I;
                if (mode == 19) T = re2 + im2 * I;
                if (mode == 20) T = re3 + im3 * I;
                if (mode == 21) T = re4 + im4 * I;
                if (mode == 22) T = re8 + im8 * I;
                if (mode == 23) T = re9 + im9 * I;
                if (mode == 24) T = re10 + im10 * I;
                if (mode == 25) T = re11 + im11 * I;
                s11ant[j] = T;
                freqs11ant[j] = freq;
                wttant[j] = 1;
                ns11ant = j + 1;
                strcpy(antfname, fname);
              }
              if (mode == 31) {  // from /data5/edges/data/CalHotLoadCableData
                T = re + im * I;
                s11rig[j] = T;
                freqs11rig[j] = freq;
                wttrig[j] = 1;
                ns11rig = j + 1;
                strcpy(rigfname, fname);
                T = re1 + im1 * I;
                s12rig[j] = T;
                freqs11rig[j] = freq;
                wttrig[j] = 1;
                ns11rig = j + 1;
                strcpy(rigfname, fname);
                T = re2 + im2 * I;
                s22rig[j] = T;
                freqs11rig[j] = freq;
                wttrig[j] = 1;
                ns11rig = j + 1;
                strcpy(rigfname, fname);
                // printf("j %d freq %f T %f
                // %s\n",j,freq,20*log10(cabs(s11rig[j])),fname);
              }
              if (mode == 33) {  // simulator 3
                s11ant[j] = re5 + im5 * I;
                freqs11ant[j] = freq;
                wttant[j] = 1;
                ns11ant = j + 1;
                strcpy(antfname, fname);
              }
              if (mode == 26) {  // cboxnew
                T = re2 + im2 * I;
                s11rig[j] = T;
                freqs11rig[j] = freq;
                wttrig[j] = 1;
                ns11rig = j + 1;
                strcpy(rigfname, fname);
                T = re3 + im3 * I;
                s12rig[j] = T;
                freqs11rig[j] = freq;
                wttrig[j] = 1;
                ns11rig = j + 1;
                strcpy(rigfname, fname);
                T = re4 + im4 * I;
                s22rig[j] = T;
                freqs11rig[j] = freq;
                wttrig[j] = 1;
                ns11rig = j + 1;
                strcpy(rigfname, fname);
                T = re1 + im1 * I;
                s11ant[j] = T;
                freqs11ant[j] = freq;
                wttant[j] = 1;
                ns11ant = j + 1;
                strcpy(antfname, fname);
              }
              if (mode == 27) {  // 3 dB for lowband2
                T = re3 + im3 * I;
                s11rig[j] = T;
                freqs11rig[j] = freq;
                wttrig[j] = 1;
                ns11rig = j + 1;
                strcpy(rigfname, fname);
                T = re4 + im4 * I;
                s12rig[j] = T;
                freqs11rig[j] = freq;
                wttrig[j] = 1;
                ns11rig = j + 1;
                strcpy(rigfname, fname);
                T = re5 + im5 * I;
                s22rig[j] = T;
                freqs11rig[j] = freq;
                wttrig[j] = 1;
                ns11rig = j + 1;
                strcpy(rigfname, fname);
                T = re1 + im1 * I;
                s11ant[j] = T;
                freqs11ant[j] = freq;
                wttant[j] = 1;
                ns11ant = j + 1;
                strcpy(antfname, fname);
              }
              if (mode == 28) {  // cbx for lowband2
                T = re10 + im10 * I;
                s11rig[j] = T;
                freqs11rig[j] = freq;
                wttrig[j] = 1;
                ns11rig = j + 1;
                strcpy(rigfname, fname);
                T = re11 + im11 * I;
                s12rig[j] = T;
                freqs11rig[j] = freq;
                wttrig[j] = 1;
                ns11rig = j + 1;
                strcpy(rigfname, fname);
                T = re12 + im12 * I;
                s22rig[j] = T;
                freqs11rig[j] = freq;
                wttrig[j] = 1;
                ns11rig = j + 1;
                strcpy(rigfname, fname);
                T = re3 + im3 * I;
                s11ant[j] = T;
                freqs11ant[j] = freq;
                wttant[j] = 1;
                ns11ant = j + 1;
                strcpy(antfname, fname);
              }
              if (mode == 29) {  // noise
                T = re + im * I;
                s11ant[j] = T;
                freqs11ant[j] = freq;
                wttant[j] = 1;
                ns11ant = j + 1;
                strcpy(antfname, fname);
              }
              j++;
            }
          }
        }
      }
      fclose(file1);
      //    printf("mode %d j %d nant %d nopen %d nhot %d namb %d ns11ant %d
      //    ns11lna %d ns11hot %d ns11amb %d ns11cab1 %d ns11cab2 %d\n",
      //          mode,j,nant,nopen,nhot,namb,ns11ant,ns11lna,ns11hot,ns11amb,ns11cab1,ns11cab2);
    }
    mode = 0;
  }
  if (tcab == 0.0) tcab = tcold;  // sets to tcold if not entered
  // wt for lnas11 needs to be decided i.e. do we want the next line
  //      for(j=0;j<ns11lna;j++) if(wtmode==0 && (freqs11lna[j] < wfstart ||
  //      freqs11lna[j] > wfstop)) wttlna[j]=0;
  for (j = 0; j < ns11lna; j++)
    s11lna[j] = s11lna[j] * cexp(-2.0 * PI * freqs11lna[j] * 1e6 * ((delaylna + delaycorr) * I)) * pow(10.0, 0.05 * (dbcorr + ldb)) + lr +
                li * I;  // corrections for adapters
  for (j = 0; j < ns11ant; j++) {
    if (opn)
      s11ant[j] = opn;
    else
      s11ant[j] = s11ant[j] * cexp(-2.0 * PI * freqs11ant[j] * 1e6 * (delaycorr + delayant) * I) * pow(10.0, 0.05 * (dbcorr + adb)) + ar + ai * I;
  }
  for (j = 0; j < ns11hot; j++) s11hot[j] = s11hot[j] * cexp(-2.0 * PI * freqs11hot[j] * 1e6 * (delaycorr)*I) * pow(10.0, 0.05 * dbcorr);
  for (j = 0; j < ns11amb; j++) s11amb[j] = s11amb[j] * cexp(-2.0 * PI * freqs11amb[j] * 1e6 * (delaycorr)*I) * pow(10.0, 0.05 * dbcorr);
  for (j = 0; j < ns11cab1; j++) s11cab1[j] = s11cab1[j] * cexp(-2.0 * PI * freqs11cab1[j] * 1e6 * (delaycorr)*I) * pow(10.0, 0.05 * dbcorr);
  for (j = 0; j < ns11cab2; j++) s11cab2[j] = s11cab2[j] * cexp(-2.0 * PI * freqs11cab2[j] * 1e6 * (delaycorr)*I) * pow(10.0, 0.05 * dbcorr);
  printf("GOT YEAR: %d %d %d %d %d\n", yr, dy, hr, mn, sc);
  if (yr)
    secs = tosecs(yr, dy, hr, mn, sc);
  else
    secs = 0;
  if (!binteg) binteg = (ghamax - ghamin) * 3600.0;
  if (binteg > 24 * 3600 || binteg < 1) binteg = 1;  // make at least 1 sec
  if (yr) sprintf(title, "%04d:%03d:%02d:%02d:%02d", yr, dy, hr, mn, sc);
  if (yrh) sprintf(hottitle, "%04d:%03d:%02d:%02d:%02d", yrh, dyh, hrh, mnh, sch);
  if (yrc) sprintf(ambtitle, "%04d:%03d:%02d:%02d:%02d", yrc, dyc, hrc, mnc, scc);
  if (yro) sprintf(opentitle, "%04d:%03d:%02d:%02d:%02d", yro, dyo, hro, mno, sco);
  if (yrs) sprintf(shorttitle, "%04d:%03d:%02d:%02d:%02d", yrs, dys, hrs, mns, scs);
  fstart = freqamb[0];
  fstop = freqamb[namb - 1];
  fcen = (fstart + fstop) / 2.0;

  if (skymode >= 0)
    bfit = beamcorr(0, bfit, antaz, secs, mcalc, fitf, aarr, bbrr, skymode, binteg, &lst, bb, &bspac, cmb, &freqref, fbstart, fbstop, site,
                    map);  // fit polynomial
  if (skymod2 >= 0)
    bfit = beamcorr2(0, bfit, antaz, secs, mcalc, fitf, aarr, bbrr, skymod2, binteg, &lst, bb, &bspac, cmb, &freqref, fbstart, fbstop,
                     site);  // fit polynomial
  if (skymod2 <= -128)
    bfit2 = beamcorr(0, bfit, antaz, secs, mcalc, fitf, aarr, bbrr, -skymod2, binteg, &lst, bb2, &bspac, cmb, &freqref, fbstart, fbstop, site,
                     map);  // fit polynomial
  if (skymode >= 0)
    sprintf(info, "lst %5.2f", lst);
  else
    sprintf(info, " ");

  if (nant && (nopen || nshort) && nhot && namb && ns11ant && ns11lna && ns11hot && ns11amb && ns11cab1 && !nocal) {  
      // continue as all data is present
      //      L = pow(10.0,-0.1*atten);  // L not used if atten needed for cable
      //      use R factor method of memo 98

      // fit s11 with series

    printf(
        "nopen %d nshort %d nhot %d namb %d ns11hot %d ns11amb %d ns11cab1 %d "
        "ns11cab2 %d freq %f %f %f %f\n",
        nopen, nshort, nhot, namb, ns11hot, ns11amb, ns11cab1, ns11cab2, freqant[0], freqant[nant - 1], freqs11ant[0], freqs11ant[ns11ant - 1]);

    // 27 or higher for antenna
    printf("nant %d ns11ant %d ns11lna %d\n", nant, ns11ant,
           ns11lna);  // period under 100 MHz
    fittp(ns11lna, freqs11lna, s11lna, wttlna, namb, freqamb, ss11lna, nfit3, fitf, mcalc, aarr, bbrr, 1, lna_poly, lnafname);
    fittp(ns11hot, freqs11hot, s11hot, wtthot, namb, freqhot, ss11hot, nfit2, fitf, mcalc, aarr, bbrr, 2, -1, hotfname);
    fittp(ns11amb, freqs11amb, s11amb, wttamb, namb, freqamb, ss11amb, nfit2, fitf, mcalc, aarr, bbrr, 3, -1, ambfname);
    if (nopen) fittp(ns11cab1, freqs11cab1, s11cab1, wttcab1, namb, freqamb, ss11cab, nfit2, fitf, mcalc, aarr, bbrr, 4, -1, openfname);
    if (nshort) fittp(ns11cab2, freqs11cab2, s11cab2, wttcab2, namb, freqamb, &ss11cab[namb], nfit2, fitf, mcalc, aarr, bbrr, 5, -1, shortfname);
    if (ns11rig) {
      fittp(ns11rig, freqs11rig, s11rig, wttrig, namb, freqamb, ss11rig, nfit2, fitf, mcalc, aarr, bbrr, 7, -1, rigfname);
      fittp(ns11rig, freqs11rig, s12rig, wttrig, namb, freqamb, ss12rig, nfit2, fitf, mcalc, aarr, bbrr, 8, -1, rigfname);
      fittp(ns11rig, freqs11rig, s22rig, wttrig, namb, freqamb, ss22rig, nfit2, fitf, mcalc, aarr, bbrr, 9, -1, rigfname);
    }

    // SGM: added for cross-checks
    FILE *s11file;
    s11file = fopen("s11_modelled.txt", "w");
    fprintf(s11file, "# freq, amb_real amb_imag hot_real hot_imag open_real open_imag short_real short_imag lna_real lna_imag rig_s11_real rig_s11_imag rig_s12_real rig_s12_imag rig_s22_real rig_s22_imag\n");
    for(i=0;i<namb;i++){
      fprintf(s11file, "%1.12e %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e\n", 
        freqhot[i],  
        creal(ss11amb[i]), cimag(ss11amb[i]),
        creal(ss11hot[i]), cimag(ss11hot[i]),
        creal(ss11cab[i]), cimag(ss11cab[i]),
        creal(ss11cab[namb + i]), cimag(ss11cab[namb + i]),
        creal(ss11lna[i]), cimag(ss11lna[i]),
        creal(ss11rig[i]), cimag(ss11rig[i]),
        creal(ss12rig[i]), cimag(ss12rig[i]),
        creal(ss22rig[i]), cimag(ss22rig[i])
      );
    }
    fclose(s11file);

    if (!nopen) {
      for (i = 0; i < namb; i++) {
        ss11cab[i] = ss11cab[i + namb];
        spopen[i] = spshort[i];
        freqopen[i] = freqshort[i];
        wtopen[i] = wtshort[i];
      }
      nopen = nshort;
      nshort = 0;  // used if only shorted cable available
    }
    // Do calibration in lab
    for (iter = 0; iter < nter; iter++) {  // needs 4 for accurate convergence
      for (i = 0; i < namb; i++) {
        if (iter == 0) {
          tcal_sca[i] = 1;
          tcal_ofs[i] = tlna0[i] = tlna1[i] = tlna2[i] = 0;
          sspopen[i] = spopen[i];
          sspopen[namb + i] = spshort[i];
          wtopen[namb + i] = wtshort[i];
          sspamb[i] = spamb[i];
          ssphot[i] = sphot[i];
          freqopen[namb + i] = freqopen[i];
        }
        // invert with scale error and offset and changed S11 (i.e. with changes
        // it will analyzed with wrong S11)
        s1 = w3pinv(ss11amb[i], ss11lna[i], sspamb[i], tload, tamb, tlna0[i], tlna1[i], tlna2[i]);
        s = w3pinv(ss11hot[i], ss11lna[i], ssphot[i], tload, tamb, tlna0[i], tlna1[i], tlna2[i]);
        if (Lh >= 0) Lhh = loss(ss11hot[i], Lh * sqrt(freqhot[i] / 100.0));
        if (Lh == -1) Lhh = lossmodel(freqhot[i], ss11hot[i],
                                      0);  // mode = 0 for cable + adapter
        if (Lh == -2) Lhh = rigloss(ss11hot[i], ss11rig[i], ss12rig[i],
                                    ss22rig[i]);  // rigid cable loss
        s2 = lossinv(s, tcold, Lhh);
        sca = (thot - tcold) / (s2 - s1);
        ofs = s1 - tcold;

        tcal_sca[i] = tcal_sca[i] * sca;
        tcal_ofs[i] += ofs;

        sspopen[i] = (spopen[i] - tamb) * tcal_sca[i] + tamb - tcal_ofs[i];
        if (nopen) sspopen[nopen + i] = (spshort[i] - tamb) * tcal_sca[i] + tamb - tcal_ofs[i];
        ssphot[i] = (sphot[i] - tamb) * tcal_sca[i] + tamb - tcal_ofs[i];
        sspamb[i] = (spamb[i] - tamb) * tcal_sca[i] + tamb - tcal_ofs[i];
        if (wtopen[i] == 0 || wthot[i] == 0 || wtamb[i] == 0 || freqamb[i] < wfstart || freqamb[i] > wfstop)
          wtcal[i] = 0;
        else
          wtcal[i] = 1;
        if (wtshort[i] == 0 || wtcal[i] == 0)
          wtcal[namb + i] = 0;
        else
          wtcal[namb + i] = 1;
        if ((wtmode % 10) == 1 && freqamb[i] >= 55.0 && freqamb[i] <= 195) wtcal[i] = wtcal[namb + i] = 1;  // allow wideband cal
        if ((wtmode % 10) == 2 && freqamb[i] >= 80.0 && freqamb[i] <= 195) wtcal[i] = wtcal[namb + i] = 1;  // allow wideband cal
        if ((wtmode % 10) == 3) {
          if (wtopen[i] == 0 || wthot[i] == 0 || wtamb[i] == 0 || freqamb[i] < 80.0 || freqamb[i] > 195.0)
            wtcal[i] = 0;
          else
            wtcal[i] = 1;
          if (wtshort[i] == 0 || wthot[i] == 0 || wtamb[i] == 0 || freqamb[i] < 80.0 || freqamb[i] > 195.0)
            wtcal[namb + i] = 0;
          else
            wtcal[namb + i] = 1;
        }
        if ((wtmode % 10) == 4) {
          if (wtopen[i] == 0 || wthot[i] == 0 || wtamb[i] == 0 || freqamb[i] < 50.0 || freqamb[i] > 100.0)
            wtcal[i] = 0;
          else
            wtcal[i] = 1;
          if (wtshort[i] == 0 || wthot[i] == 0 || wtamb[i] == 0 || freqamb[i] < 50.0 || freqamb[i] > 100.0)
            wtcal[namb + i] = 0;
          else
            wtcal[namb + i] = 1;
        }
        if (!(i % (nopen / 10)) && iter >= nter - 1)
          printf(
              "iter %d freq %7.2f sca %7.2f s1 %7.2f s %7.2f s2 %7.2f tcold "
              "%7.2f thot %7.2f ofs %5.2f sca %5.2f tcal_ofs %5.2f tcal_sca "
              "%5.2f tlnau %5.2f\n",
              iter, freqhot[i], sca, s1, s, s2, tcold, thot, ofs, sca, tcal_ofs[i], tcal_sca[i], tlna0[i]);
      }
      //    for(i=0;i<namb;i++) Tll[i]=Tll2[i]=s11lna[i]; // test of fitting
      wavefit(nopen, nshort, wfit, sspopen, wtcal, freqopen, ss11cab, ss11lna, tcab, tload, tamb, tlna0, tlna1, tlna2, dataout, mcalc, fitf, aarr,
              bbrr, &delayln);
    }
    printf("Got best delay: %g\n", delayln);

    plotfspec(nopen + nshort, freqopen, sspopen, dataout, wtcal, dscale, 0, 3, openfile, shortfile, info);  // check on fit to open cable
    plotnwave(namb, freqopen, tlna0, tlna1, tlna2, delayln, wfit);
    printf("rms %f K wfit %d nopen %d\n", rmscalc(nopen + nshort, sspopen, dataout, wtcal), wfit, nopen + nshort);

    for (i = 0; i < namb; i++) {
      data[i] = thot;
      s2 = w3pinv(ss11hot[i], ss11lna[i], ssphot[i], tload, tamb, tlna0[i], tlna1[i], tlna2[i]);
      if (Lh < 0)
        Lhh = lossmodel(freqhot[i], ss11hot[i], 0);
      else
        Lhh = loss(ss11hot[i], Lh * sqrt(freqhot[i] / 100.0));
      dataout[i] = lossinv(s2, tcold, Lhh);
    }  // note tcold
    printf("rms %f K wfit %d nhot %d\n", rmscalc(namb, data, dataout, wtcal), wfit, nhot);
    //   plotfspec(namb,freqhot, sphot, dataout,
    //   wtcal,dscale,0,1,title,datfile,info);
    for (i = 0; i < namb; i++) {
      data[i] = tcold;
      dataout[i] = w3pinv(ss11amb[i], ss11lna[i], sspamb[i], tload, tamb, tlna0[i], tlna1[i], tlna2[i]);
    }
    printf("rms %f K wfit %d namb %d\n", rmscalc(namb, data, dataout, wtcal), wfit, namb);
    //   plotfspec(namb,freqamb, data, dataout,
    //   wtcal,dscale,0,2,title,datfile,info);

    // 7 terms to smooth cal corrections
    for (j = 0; j < cfit; j++) {
      for (i = 0; i < namb; i++) {
        fitf[i + j * namb] = pow((freqamb[i] / fcen),
                                 0.5 * j);  // poly is best cfit 6 - removed Fourier
                                            // series - add 0.5 may be better
      }
    }

    if (cfit > 0) {
      polyfitr(cfit, namb, tcal_sca, mcalc, wtcal, tcal_scasm, fitf, aarr, bbrr);
      printf("rms_tcal_sca %f dB cfit %d\n", 4.0 * rmscalc(namb, tcal_sca, tcal_scasm, wtcal),
             cfit);  // scale to get in dB approx
      polyfitr(cfit, namb, tcal_ofs, mcalc, wtcal, tcal_ofssm, fitf, aarr, bbrr);
      printf("rms_tcal_ofs %f K cfit %d\n", rmscalc(namb, tcal_ofs, tcal_ofssm, wtcal), cfit);
      //   for(i=0;i<namb;i++) printf("freq %f tcal_ofs %f tcal_ofssm %f diff %f
      //   wtcal
      //   %f\n",freqamb[i],tcal_ofs[i],tcal_ofssm[i],tcal_ofs[i]-tcal_ofssm[i],wtcal[i]);
    } else
      for (i = 0; i < namb; i++) {
        tcal_scasm[i] = tcal_sca[i];
        tcal_ofssm[i] = tcal_ofs[i];
      }  // no fitting
    outcal(namb, freqamb, ss11lna, tcal_scasm, tcal_ofssm, tlna0, tlna1, tlna2, wtcal, "cal_data");
    cal = 1;

    FILE *LOSS_FILE;
    LOSS_FILE = fopen("hot_load_loss.txt", "w");

    fprintf(LOSS_FILE, "# freq, loss\n");
  
    for (i = 0; i < namb; i++) {
      data[i] = thot;
      ssphot[i] = (sphot[i] - tamb) * tcal_scasm[i] + tamb - tcal_ofssm[i];
      s2 = w3pinv(ss11hot[i], ss11lna[i], ssphot[i], tload, tamb, tlna0[i], tlna1[i], tlna2[i]);
      if (Lh >= 0) Lhh = loss(ss11hot[i], Lh * sqrt(freqhot[i] / 100.0));
      if (Lh == -1) Lhh = lossmodel(freqhot[i], ss11hot[i],
                                    0);  // mode = 0 for cable + adapter
      if (Lh == -2) Lhh = rigloss(ss11hot[i], ss11rig[i], ss12rig[i],
                                  ss22rig[i]);  // rigid cable loss
      fprintf(LOSS_FILE, "%1.12e %1.12e\n", freqhot[i], Lhh);
      dataout[i] = lossinv(s2, tcold, Lhh);     // note tcold
    }
    fclose(LOSS_FILE);
    plotfspec(namb, freqhot, dataout, data, wtcal, dscale, 0, 1, hottitle, hotfile, info);

    for (i = 0; i < namb; i++) {
      data[i] = tcold;
      sspamb[i] = (spamb[i] - tamb) * tcal_scasm[i] + tamb - tcal_ofssm[i];
      dataout[i] = w3pinv(ss11amb[i], ss11lna[i], sspamb[i], tload, tamb, tlna0[i], tlna1[i], tlna2[i]);
    }
    plotfspec(namb, freqamb, dataout, data, wtcal, dscale, 0, 2, ambtitle, ambfile, info);

  } else if (!nant) {
    printf("missing data file\n");
    return 0;
  }
  if (cal == 0 || nant != namb) {
    readcal(&ncal, freqcal, ss11lna, tcal_scasm, tcal_ofssm, tlna0, tlna1, tlna2, wtcal);
    if (ncal == 0) return 0;
    for (i = 0; i < ncal; i++) tmp[i] = ss11lna[i] * cexp(-2.0 * PI * freqcal[i] * 1e6 * (delaylna * I)) * pow(10.0, 0.05 * ldb) + lr + li * I;
    fittp(ncal, freqcal, tmp, wtcal, nant, freqant, tmpsm, nfit1, fitf, mcalc, aarr, bbrr, 6, -1, fname);
    for (i = 0; i < nant; i++) ss11lna[i] = tmpsm[i];
    for (i = 0; i < ncal; i++) tmp[i] = tcal_scasm[i] + tcal_ofssm[i] * I;
    fittp(ncal, freqcal, tmp, wtcal, nant, freqant, tmpsm, nfit1, fitf, mcalc, aarr, bbrr, 6, -1, fname);
    for (i = 0; i < nant; i++) {
      tcal_scasm[i] = creal(tmpsm[i]);
      tcal_ofssm[i] = cimag(tmpsm[i]);
    }
    for (i = 0; i < ncal; i++) tmp[i] = tlna1[i] + tlna2[i] * I;
    fittp(ncal, freqcal, tmp, wtcal, nant, freqant, tmpsm, nfit1, fitf, mcalc, aarr, bbrr, 6, -1, fname);
    for (i = 0; i < nant; i++) {
      tlna1[i] = creal(tmpsm[i]);
      tlna2[i] = cimag(tmpsm[i]);
    }
    for (i = 0; i < ncal; i++) tmp[i] = tlna0[i];
    fittp(ncal, freqcal, tmp, wtcal, nant, freqant, tmpsm, nfit1, fitf, mcalc, aarr, bbrr, 6, -1, fname);
    for (i = 0; i < nant; i++) tlna0[i] = creal(tmpsm[i]);
    for (i = 0; i < nant; i++)
      if (freqant[i] < freqcal[0] || freqant[i] > freqcal[ncal - 1]) wtant[i] = 0;
  }
  if (wtmode > 10)
    for (i = 0; i < nant; i++)
      if (freqs11ant[i] < wfstart - 1 || freqs11ant[i] > wfstop + 1) wttant[i] = 0;  // useful for lower order fit
  kk = 0;
  for (j = 0; j < ns11ant; j++)
    if (wttant[j]) kk++;
  if (nfit4 >= kk) nfit4 = kk - 1;
  fittp(ns11ant, freqs11ant, s11ant, wttant, nant, freqant, ss11ant, nfit4, fitf, mcalc, aarr, bbrr, 0, -1, antfname);
  if (lmode == 2) {
    fittp(ns11rig, freqs11rig, s11rig, wttrig, nant, freqant, ss11rig, nfit2, fitf, mcalc, aarr, bbrr, 7, -1, rigfname);
    fittp(ns11rig, freqs11rig, s12rig, wttrig, nant, freqant, ss12rig, nfit2, fitf, mcalc, aarr, bbrr, 8, -1, rigfname);
    fittp(ns11rig, freqs11rig, s22rig, wttrig, nant, freqant, ss22rig, nfit2, fitf, mcalc, aarr, bbrr, 9, -1, rigfname);
  }


  if ((freqant[nant - 1] + freqant[0]) / 2.0 < 100.0)
    freqr = 75.0;
  else
    freqr = 150.0;
  if (freqref > freqr * 1.5) low = 1;
  
  FILE *lossfile;
  FILE *beamcorrfile;
  double refsky;
  
  lossfile = fopen("loss.txt", "w");
  beamcorrfile = fopen("beamcorr.txt", "w");
  
  fprintf(lossfile, "# freq, loss, tloss [(1 - loss)*Tant]\n");
  fprintf(beamcorrfile, "# freq, skymodel, refskymodel, beamcorr\n");

  printf("WTANT JUST BEFORE LOSS/BEAM: %e\n", wtant[0]);

  for (i = 0; i < nant; i++) {
    sspant[i] = (spant[i] - tamb) * tcal_scasm[i] + tamb - tcal_ofssm[i];
    //   if(freqant[i] < wfstart || freqant[i] > wfstop || wtcal[i] == 0)
    //   wtant[i] = 0;
    if (freqant[i] < wfstart || freqant[i] > wfstop) wtant[i] = 0;
    dataout[i] = w3pinv(ss11ant[i], ss11lna[i], sspant[i], tload, tamb, tlna0[i], tlna1[i], tlna2[i]);

    if (lmode < 0)
      La = loss(ss11ant[i],
                atten + aloss * sqrt(freqant[i] / 100.0));  // aloss > 0 for simple atten or cross-check
    else {                                                  // 6 = no antenna or ground plane loss
      if (lmode == 2)
        La = rigloss(ss11ant[i], ss11rig[i], ss12rig[i],
                     ss22rig[i]);  // hot load - can only use if calibration
                                   // being done otherwise ss11rig=0
      else
        La = lossmodel(freqant[i], ss11ant[i], lmode) * antloss(freqant[i], mcalc, fitf, aarr, bbrr, aloss,
                                                                lmode);  // lmode 0=hot load 1=balun+bend 5or6=balun+bend
                                                                         // lowband
    }
    if (wtant[i] && (La < 0.1 || La > 1.0)) printf("check loss freq %f La %f\n", freqant[i], La);

    ltemp[i] = (1.0 - La) * tant;

    fprintf(lossfile, "%f %e %e\n", freqant[i], La, ltemp[i]);

    dataout[i] = lossinv(dataout[i], tant, La);

    if (!(i % (nant / 10)))
      printf(
          "freqant %7.2f tlnau %7.2f tlnacc %7.2f tlnacs %7.2f tlnac %7.2f "
          "phase %7.2f tsky %7.2f md %7.2f loss %5.2f\n",
          freqant[i], tlna0[i], tlna1[i], tlna2[i], sqrt(tlna1[i] * tlna1[i] + tlna2[i] * tlna2[i]), atan2(tlna2[i], tlna1[i]) * 180.0 / PI,
          dataout[i], md[i], (1.0 - La) * tant);
    //        printf("loss %f dB La %f freq %f dataout %f md %f beamcr
    //        %f\n",loss,La,freqant[i],dataout[i],md[i],beamcr); beamcorr
    //        assumes a spectral index of -2.5
    if (test & 16 && skymod2 <= -128)
      dataout[i] =
          beamcorr(freqant[i] * (1.0 + low), bfit2, antaz, secs, mcalc, fitf, aarr, bbrr, 0, 0, &lst, bb2, &bspac, cmb, &freqref, 0, 0, site, map) *
          pow(1.0 / (1.0 + low), -2.5);  // substitute  beamcorrection
    
    if ((skymode >= 0) || (skymod2 >= 0)){
      if (skymode >= 0) {
        skymodel[i] =
            beamcorr(freqant[i] * (1.0 + low), bfit, antaz, secs, mcalc, fitf, aarr, bbrr, 0, 0, &lst, bb, &bspac, cmb, &freqref, 0, 0, site, map) *
            pow(1.0 / (1.0 + low), -2.5);
        refsky = 
            beamcorr(freqref, bfit, antaz, secs, mcalc, fitf, aarr, bbrr, 0, 0, &lst, bb, &bspac, cmb, &freqref, 0, 0, site, map) *
            pow(freqant[i] / freqref, -2.5);
      }
      if (skymod2 >= 0) {
        skymodel[i] =
            beamcorr2(freqant[i] * (1.0 + low), bfit, antaz, secs, mcalc, fitf, aarr, bbrr, 0, 0, &lst, bb, &bspac, cmb, &freqref, 0, 0, site) *
            pow(1.0 / (1.0 + low), -2.5);
        refsky = 
            beamcorr2(freqref, bfit, antaz, secs, mcalc, fitf, aarr, bbrr, 0, 0, &lst, bb, &bspac, cmb, &freqref, 0, 0, site) *
            pow(freqant[i] / freqref, -2.5) ;
      }

      dataout[i] = dataout[i] * refsky / skymodel[i];

      fprintf(beamcorrfile, "%f %e %e %e\n", freqant[i], skymodel[i], refsky, refsky / skymodel[i]);
    }

    if (test & 1 && (skymode >= 0 || skymod2 >= 0) && (skymode >= 256 || skymod2 >= 256)) dataout[i] = skymodel[i];  // substitute  beamcorrection
    if (test == 2) dataout[i] = 300 * pow(freqant[i] / 150.0,
                                          specin);  // test rms 14 mK at 2.55 1 mK at 2.51
    if (test == 4)
      dataout[i] = 300 * pow(freqant[i] / 150.0,
                             -2.55 + 0.01 * log(freqant[i] / 150.0));  // test rms 14 mK at 2.55 1
                                                                       // mK at 2.51 gamma 0.01
    if (test == 6) {
      dataout[i] = 300 * pow(freqant[i] / 150.0,
                             -2.55 + 0.0 * log(freqant[i] / 150.0));  // gamma 0.0
      dataout[i] += dataout[i] * (-1e-3 * pow(freqant[i] / 150.0, -2.0)) + 1e3 * 1e-3 * pow(freqant[i] / 150.0, -2.0);
    }  // 0.1% ion at 1000K
    if (test & 8) {
      if (tau > 0)
        dataout[i] +=
            -eoramp *
            (1 -
             exp(-tau * exp(-(freqant[i] - eorcen) * (freqant[i] - eorcen) * (-log(-log((1 + exp(-tau)) / 2) / tau)) / (eorwid * eorwid * 0.25)))) /
            (1 - exp(-tau));
      else
        dataout[i] += -eoramp * exp(-0.69 * (freqant[i] - eorcen) * (freqant[i] - eorcen) / ((eorwid * eorwid + 1e-6) * 0.25));
      dataout[i] += dataout[i] * ion * (-1e-3 * pow(freqant[i] / 150.0, -2.0)) + ion * 1e3 * 1e-3 * pow(freqant[i] / 150.0, -2.0);
    }  // 1% ion at 1000K
  }
  fclose(lossfile);
  fclose(beamcorrfile);

  if (sim) {
    for (i = 0; i < nant; i++) {
      if (test & 1 && (skymode >= 0 || skymod2 >= 0) && (skymode >= 256 || skymod2 >= 256)) {
        tant2 = skymodel[i] + skymodel[i] * ion * (-1e-3 * pow(freqant[i] / 150.0, -2.0)) + ion * 1e3 * 1e-3 * pow(freqant[i] / 150.0, -2.0);
        if (tau > 0)
          tant2 +=
              -eoramp *
              (1 -
               exp(-tau * exp(-(freqant[i] - eorcen) * (freqant[i] - eorcen) * (-log(-log((1 + exp(-tau)) / 2) / tau)) / (eorwid * eorwid * 0.25)))) /
              (1 - exp(-tau));
        else
          tant2 += -eoramp * exp(-0.69 * (freqant[i] - eorcen) * (freqant[i] - eorcen) / ((eorwid * eorwid + 1e-6) * 0.25));
      } else
        tant2 = sim * 300.0 * pow(freqant[i] / 150.0, -2.5);
      if (sim % 10 == 9) tant2 = dataout[i];
      La = 1.0 - ltemp[i] / tant;
      if (sim / 10) La = 1;  // no loss correction
      if (sim / 10 == 2) La = lossmodel(freqant[i], ss11ant[i], lmode);
      tant2 = La * tant2 + (1.0 - La) * tant;
      data[i] = (w3p(ss11ant[i], ss11lna[i], tant2, tload, tamb, tlna0[i], tlna1[i], tlna2[i]) - tamb + tcal_ofssm[i]) / tcal_scasm[i] + tamb +
                noise * gauss();
    }
    outsim(nant, freqant, data, yr, dy, hr, mn, sc);
  }

  if (ltemp[0]) plotfspec(nant, freqant, ltemp, ltemp, wtant, dscale, 0, 5, title, datfile, info);
  plotfspec(nant, freqant, md, dataout, wtant, dscale, 0, 0, title, datfile, info);

  rms = rms1 = rmscalc(nant, dataout, md, wtant);
  err1 = 0;
  i1 = i2 = i3 = i4 = i5 = i6 = i7 = -1;  // none selected
  for (kk = 0; kk < 2; kk++) {
    if (kk == 0)
      nfit = 2;
    else
      nfit = abs(mfit);
    for (i = 0; i < nfit; i++) {
      for (j = 0; j < nant; j++) {
        k = j + i * nant;
        freq = freqant[j];
        fcen = (freqant[nant - 1] + freqant[0]) / 2.0;
        if (eorcen == 0) {
          if (i == 0) fitf[k] = 1;  // constant
          if (i > 0 && (i % 2) == 1) {
            fitf[k] = w3pinv(ss11ant[j], ss11lna[j], sspant[j], tload, tamb, tlna0[j], tlna1[j],
                             tlna2[j])  // s11ant mag
                      - w3pinv(ss11ant[j] * pow(10.0, -0.05 * 0.01 * pow(freq / fcen, (i / 2))), ss11lna[j], sspant[j], tload, tamb, tlna0[j],
                               tlna1[j], tlna2[j]);
          }

          if (i > 0 && (i % 2) == 0) {
            fitf[k] = w3pinv(ss11ant[j], ss11lna[j], sspant[j], tload, tamb, tlna0[j], tlna1[j],
                             tlna2[j])  // s11ant phase
                      - w3pinv(ss11ant[j] * cexp(-(PI / 180.0) * pow(freq / fcen, (i / 2) - 1) * I), ss11lna[j], sspant[j], tload, tamb, tlna0[j],
                               tlna1[j], tlna2[j]);
          }
        } else {
          if (mdd == 1) {
            i1 = 1;
            i3 = 2;
            i4 = 3;
          }  // 4 terms sca,sp_ind,ion_abs,ion_em
          if (mdd == 2) {
            i3 = 1;
            i4 = 2;
          }  // 3 terms sca,ion_abs,ion_em
          if (mdd == 4) {
            i1 = 1;
            i2 = 2;
            i3 = 3;
            i4 = 4;
            i5 = 5;
          }  // 5 (6) terms sca,sp_ind,gamma,ion_abs,ion_em,(cons)
          if (mdd == 5) {
            i3 = 1;
            i6 = 2;
          }                          // 3 terms sca,ion_abs,sun
          if (mdd == 6) { i6 = 1; }  // 2 terms sca,sun
          if (mdd == 7) {
            i3 = 1;
            i4 = 2;
            i6 = 3;
          }  // 4 terms sca,ion_abs,ion_em,sun
          if (mdd == 8) {
            i1 = 1;
            i2 = 2;
            i7 = 3;
            i3 = 4;
            i4 = 5;
            i5 = 6;
          }  // 5(6) terms sca,sp_ind,gamma,curve,sca,ion_abs,ion_em
          if (i == 0) {
            if (mfit > 0)
              fitf[k] = pow(freq / freqr, specin);
            else
              fitf[k] = 1;
          }                                                                      // scale mfit=-1 useful for test
          if (i == i1) fitf[k] = log(freq / freqr) * pow(freq / freqr, specin);  // adjustment of spec. index
          if (i == i2)
            fitf[k] = log(freq / freqr) * log(freq / freqr) * pow(freq / freqr, specin);  // slope of index "gamma"
                                                                                          //       if(i==i2) fitf[k] =
          //       log(freq/freqr)*((freq-freqr)/freqr)*pow(freq/freqr,specin);
          //       linear case almost same
          if (i == i7)
            fitf[k] = pow(log(freq / 150.0), 3.0) * pow(freq / 150.0, specin);  // next term
                                                                                //       if(i==4) fitf[k] =
          //       pow(log(freq/150.0),4.0)*pow(freq/150.0,specin);   // next
          //       term
          if (i == i3)
            fitf[k] = pow(freq / freqr,
                          -2.0 + specin);  // ion absorption should be negative
                                           // and multiplicative
          if (i == i4) fitf[k] = pow(freq / freqr,
                                     -2.0);  // ion emission should be positive and additive
          if (i == i5) fitf[k] = 1;          // constant
          //       ss0=2.8; ss1=2.15;  var0=0.1; var1=0.01;
          if (i == i6)
            fitf[k] = pow(freq / freqr,
                          sunind);  // Sun spe_ind 0 (quiet to -2(activ))   -
                                    // see Kraus  0.5 best
          //       if(i==10) fitf[k] = 1.0/(1-ss11ant[j]*conj(ss11ant[j])); //
          //       some added instrumental functions will be needed for ant. s11
          //       structure
          if (!mdd)
            fitf[k] = pow(freq / freqr,
                          specin + i);  // best according to tests 5 + 1 terms
                                        // 20 MHz full-width err=2 40 MHz err=6
                                        //       if(i==11) fitf[k] = 1.0 -
          //       loss(ss11ant[j],0,0,0,0,freqant[i],dlyrob,lmode); // some
          //       added instrumental functions will be needed for ant. s11
          //       structure if((i==nfit-1 && nfit > 1 && eorwid > 0.0) || (i==1
          //       && kk==0 && eorwid > 0.0)) fitf[k] =
          //       exp(-0.69*(freq-eorcen)*(freq-eorcen)/(eorwid*eorwid));  //
          //       eor
        }
      }
    }
    for (i = 0; i < 40; i++) bbrr[i] = 0;  // zero previous soln.
    if (cmb & 1)
      for (i = 0; i < nant; i++) dataout[i] += -2.725;  // subtract cmb for fit
    polyfitrc(nfit, nant, dataout, mcalc, wtant, dataout2, fitf, aarr, bbrr, cons);
    if (cmb & 1)
      for (i = 0; i < nant; i++) {
        dataout[i] += 2.725;
        dataout2[i] += 2.725;
      }  // add back cmb
    if (kk == 0) err1 = sqrt(aarr[nfit - 1 + (nfit - 1) * nfit]);
  }

  FILE *ants11file;
  ants11file = fopen("modeled_antenna_s11.txt", "w");
  fprintf(ants11file, "# freq, re(s11), im(s11)\n");
  for (i=0; i < nant; i++){
    fprintf(ants11file, "%f %e %e\n",
      freqant[i], 
      creal(ss11ant[i]), cimag(ss11ant[i])
    );
  }
  fclose(ants11file);

  if (smooth) {
    FILE *smoothinputfile;
    smoothinputfile = fopen("second_smooth_input.txt", "w");
    fprintf(smoothinputfile, "# freq, data, model, wt\n");
    for (i=0; i < nant; i++){
      fprintf(smoothinputfile, "%f %e %e %e\n",
        freqant[i], 
        dataout[i], dataout2[i], wtant[i]
      );
    }
    fclose(smoothinputfile);
    dsmooth(nant, dataout, dataout2, wtant, mcalc, fabs(smooth));
    if (smooth < 0) {
      j = 0;
      kk = fabs(smooth);
      for (i = 0; i < nant; i += kk) {
        wtemp[j] = 0;
        for (k = j * kk - kk / 2; k <= j * kk + kk / 2; k++)
          if (k >= 0 && k < nant && wtant[k]) wtemp[j] = 1;

        dataout[j] = dataout[i];
        freqant[j] = freqant[i];
        mcalc[j] = mcalc[i];
        dataout2[j] = dataout2[i];
        wtant[j] = wtemp[j];
        j++;
      }
      nant = j;
    }
  } else
    for (i = 0; i < nant; i++) mcalc[i] = dataout[i];
  
  FILE *smoothoutputfile;
  smoothoutputfile = fopen("second_smooth_output.txt", "w");
  fprintf(smoothoutputfile, "# freq, data, model, wt\n");
  for (i=0; i < nant; i++){
    fprintf(smoothoutputfile, "%e %e %e %e\n",
      freqant[i], 
      mcalc[i], dataout2[i], wtant[i]
    );
  }
  fclose(smoothoutputfile);

  rms1 = rms2 = rms3 = rmscalc(nant, mcalc, dataout2, wtant);
  if (smooth) rms3 = rms3 * sqrt(fabs(smooth));
  if (eorwid > 0.0 && cmb & 32) {
    for (i = 0; i < nant; i++) wtemp[i] = wtant[i] * exp(-0.69 * (freqant[i] - eorcen) * (freqant[i] - eorcen) / (eorwid * eorwid));
    rms2 = rmscalc(nant, mcalc, dataout2, wtemp);
    if (smooth)
      rms2 = rms2 * sqrt(fabs(smooth));  // to correct correlated noise since
                                         // smooth does not reduce number of freqs
  }
  t150 = 0;
  for (i = 0; i < nfit; i++)
    if (!mdd || (i != i1 && i != i2)) t150 += bbrr[i];
  i = nfit - 1;
  if (eorcen > 100.0)
    days = aarr[i + i * nfit] * 1e2 * 300.0 * 300.0 / (20e-3 * 20e-3 * (freqant[1] - freqant[0]) * 1e6 * 3600.0 * 8.0 * 0.3 * 0.4 * 0.5);
  else
    days = aarr[i + i * nfit] * 1700.0 * 1700.0 * 1e2 /
           (100e-3 * 100e-3 * (freqant[1] - freqant[0]) * 1e6 * 3600.0 * 8.0 * 0.3 * 0.4 * 0.5);  // 8hrs,3pos,40%eff,loadnoise
  if (eorwid > 0.0)
    printf(
        "Summary mfit %d mdd %d rms %4.0f K rms1 %6.1f mK value %8.1f err "
        "%5.1f snr %3.0f days %3.0f gha %6.2f T%d %5.1f sqrtcov %5.1f\n",
        mfit, mdd, rms, 1e3 * rms1, 1e3 * bbrr[i], sqrt(aarr[i + i * nfit]) / err1, bbrr[i] / (sqrt(aarr[i + i * nfit]) * rms2), days, gha,
        (int)freqr, t150, sqrt(aarr[i + i * nfit]));
  else
    printf(
        "Summary mfit %d mdd %d rms %4.0f K rms1 %6.1f mK gha %6.2f T%d "
        "%5.1f\n",
        mfit, mdd, rms, 1e3 * rms1, gha, (int)freqr, t150);
  gamma = ion_abs = ion_em = ssun = spe_inderr = gammaerr = 0;
  spe_ind = specin;
  if (i1 > 0 && i1 < nfit) {
    spe_ind = specin + bbrr[i1] / bbrr[0];
    spe_inderr = rms3 * sqrt(aarr[i1 + i1 * nfit]) / bbrr[0];
  }
  if (i3 > 0 && i3 < nfit) ion_abs = -bbrr[i3] / bbrr[0];
  if (i4 > 0 && i4 < nfit)
    if (ion_abs) ion_em = bbrr[i4] / ion_abs;
  if (i2 > 0 && i2 < nfit) {
    gamma = bbrr[i2] / bbrr[0];
    gammaerr = rms3 * sqrt(fabs(aarr[i2 + i2 * nfit])) / bbrr[0];
  }
  // note errors are not realistic unless random noise
  if (i6 > 0 && i6 < nfit) ssun = bbrr[i6];
  if ((i1 > 0 && i1 < nfit) || (i3 > 0 && i3 < nfit))
    printf(
        "Spectral_index %6.3f +/- %5.3f gamma %6.3f +/- %5.3f T%d %5.1f "
        "ion_opac %5.3f ion_em %5.0f K rms %6.3f K mfit %d mdd %d date "
        "%04d:%03d:%02d:%02d:%02d gha %8.3f\n",
        spe_ind, spe_inderr, gamma, gammaerr, (int)freqr, t150, ion_abs, ion_em, rms1, mfit, mdd, yr, dy, hr, mn, sc, gha);
  plotfspec(nant, freqant, mcalc, dataout2, wtant, dscale, 0, 99, title, datfile, info);
  outfspec(nant, freqant, mcalc, dataout2, wtant, 0, 0, nline, title);
  return 0;
}

void wavefit(int nopen, int nshort, int wfit, double spopen[], double wtopen[], double freqopen[], complex double s11cab[], complex double s11lna[],
             double tcab, double tload, double tamb, double tlna0[], double tlna1[], double tlna2[], double dataout[], double mcalc[], double fitf[],
             long double aarr[], double bbrr[], double *delay) {
  int i, j, k, nfit, m, mopen, iter;
  double freq, tlnau, tlnacc, tlnacs, cc, cs, dly, min, bdly, rms, dl, dd;
  static double dataout2[NDATA];
  nfit = 3 * wfit;              // 9 is too many doesn't help fit
  tlnau = tlnacc = tlnacs = 0;  // dly=0;  // LNA delay
  mopen = nopen + nshort;
  rms = 0;
  min = 1e99;
  bdly = 0;
  
  for (iter = 0; iter < 2; iter++) {
    if (iter == 0)
      dd = 10e-9;
    else
      dd = 0;
    for (dl = 0; dl <= dd; dl += 1e-9) {  // not clear a delay search is needed
      if (iter == 0)
        dly = -dl * 1e6;
      else {
        dly = -bdly * 1e6;
        *delay = bdly;
      }
      for (j = 0; j < nfit; j++) {
        for (i = 0; i < mopen; i++) {
          m = i % nopen;
          freq = freqopen[m] / freqopen[nopen / 2];
          if (j % 3 == 0) {
            tlnau = pow(freq, j / 3);
            tlnacc = 0.0;
            tlnacs = 0;
          }
          if (j % 3 == 1) {
            tlnau = 0.0;
            tlnacc = pow(freq, j / 3);
            tlnacs = 0;
          }
          if (j % 3 == 2) {
            tlnau = 0.0;
            tlnacc = 0.0;
            tlnacs = pow(freq, j / 3);
          }
          cc = tlnacc * cos(2.0 * PI * freqopen[m] * dly) + tlnacs * sin(2.0 * PI * freqopen[m] * dly);
          cs = tlnacs * cos(2.0 * PI * freqopen[m] * dly) - tlnacc * sin(2.0 * PI * freqopen[m] * dly);
          fitf[i + j * mopen] =
              w3p(s11cab[i], s11lna[m], tcab, tload, tamb, tlnau, cc, cs) - w3p(s11cab[i], s11lna[m], tcab, tload, tamb, 0, 0,
                                                                                0);  // need to subract value for tlnau,tlnacc,tlnacs all zero
          dataout[i] = spopen[i] - w3p(s11cab[i], s11lna[m], tcab, tload, tamb, 0, 0,
                                       0);  // need to subtract cable noise
        }
      }
      polyfitr(nfit, mopen, dataout, mcalc, wtopen, dataout2, fitf, aarr, bbrr);
      rms = rmscalc(mopen, dataout, dataout2, wtopen);
      //  printf("Sp rms %f\n",rms);
      if (rms < min) {
        min = rms;
        bdly = dl;
      }
    }
  }
  for (i = 0; i < mopen; i++) {
    m = i % nopen;
    freq = freqopen[m] / freqopen[nopen / 2];
    tlnau = tlnacc = tlnacs = 0;
    for (k = 0; k < nfit / 3; k++) {
      tlnau += bbrr[3 * k] * pow(freq, k);
      tlnacc += bbrr[3 * k + 1] * pow(freq, k);
      tlnacs += bbrr[3 * k + 2] * pow(freq, k);
    }
    cc = tlnacc * cos(2.0 * PI * freqopen[m] * dly) + tlnacs * sin(2.0 * PI * freqopen[m] * dly);
    cs = tlnacs * cos(2.0 * PI * freqopen[m] * dly) - tlnacc * sin(2.0 * PI * freqopen[m] * dly);
    dataout[i] = w3p(s11cab[i], s11lna[m], tcab, tload, tamb, tlnau, cc,
                     cs);  // treat like an antenna
    tlna0[i] = tlnau;
    tlna1[i] = cc;
    tlna2[i] = cs;
    //   if(!(i%10))printf("freq %f data %f S11 %f
    //   %f\n",freqopen[i],dataout[i],creal(s11cab[i]),cimag(s11cab[i]));
  }
}

void fittp(int np, double freqq2[], complex double s11[], double wtt[], int npout, double freqq1[], complex double Taa[], int nfit, double fitf[],
           double mcalc[], long double aarr[], double bbrr[], int mode, int model, char *fname) {
  /*
      2022-24-03 Steven Murray: Added the "model" parameter, which specifies either

        < 0: continue to do what happened in Alan's original code (i.e. use Fourier
             when nfit > 16 and polynomial otherwise)
        0:   use Fourier
        1:   use Polynomial

      This was added so that we could switch back to previous behaviour of Alan's code
      in which 11 terms would also use a Fourier series.
  */
  int i, j, k, i1, i2;
  double data[NDATA], dataout[NDATA], dbdiff[NDATA], pdiff[NDATA], sum, a, dfreq, rmsdb, rmsphase, f0, fcen, max, del, delay;
  complex double Taaa[NDATA], sumc;
  delay = 0;
  if (mode < 6) {  // find best delay
    max = -1e99;
    for (del = -1e-9; del <= 100e-9; del += 1e-10) {
      sumc = 0;
      for (i = 0; i < np; i++) { sumc += wtt[i] * s11[i] * cexp(2.0 * PI * freqq2[i] * 1e6 * del * I); }
      if (cabs(sumc) > max) {
        max = cabs(sumc);
        delay = del;
      }
    }
    printf("Sp delay for mode %d: %e\n", mode, delay);
  }

  // SGM: modified the next two lines to start with i1 = -1 -- otherwise you always
  //      get at least i1=1 (not zero), which means you never fit down to the lowest
  //      frequency.
  // i1 = 0;
  // for(i=0;i<np && i1==0;i++) if(wtt[i]) i1=i;
  i1 = -1;
  for(i=0;i<np && i1==-1;i++) if(wtt[i]) i1=i;
  // SGM: END ADDITION

  i2 = 0;
  for (i = np - 1; i >= 0 && i2 == 0; i--)
    if (wtt[i]) i2 = i;
  f0 = freqq2[i1];
  dfreq = (freqq2[i2] - f0) * 1.5;  // was 2 might be better also wt = 1/cabs(T)
  if (dfreq == 0.0) {
    printf("no data in fittp\n");
    return;
  }
  fcen = (freqq2[i2] + f0) * 0.5;  // fixed 25sep15 now O.K. to nfit4 = 10
  

  
  for (j = 0; j < nfit; j++) {
    for (i = 0; i < np; i++) {
      if ((nfit > 16)  | (model == 0)) {  // changed 13dec17 to 14 7aug18 to 16 20aug18
        if (j == 0) fitf[i + j * np] = 1;
        if (j % 2) fitf[i + j * np] = cos((freqq2[i] - f0) * 2.0 * PI * (j / 2 + 1) / dfreq);
        if (j % 2 == 0 && j > 0) fitf[i + j * np] = sin((freqq2[i] - f0) * 2.0 * PI * (j / 2) / dfreq);
      } else if ( (model == 1) | (model < 0)) {  //
        // log  helps a little  pow(log10(freqq2[i]/fcen),j);
        fitf[i + j * np] = pow(log10(freqq2[i] / fcen), j);
      } else {
        printf("model must be -1, 0 or 1\n");
        return;
      }
    }
  }
  for (k = 0; k < 2; k++) {
    for (i = 0; i < np; i++) {
      if (k == 0)
        data[i] = creal(s11[i] * cexp(2.0 * PI * freqq2[i] * 1e6 * delay * I));
      else
        data[i] = cimag(s11[i] * cexp(2.0 * PI * freqq2[i] * 1e6 * delay * I));
    }
    // if(mode==2) for(i=0;i<np;i++) printf("k %d i %d data %f s11 %f %f delay
    // %f freqq2 %f\n",k,i,data[i],creal(s11[i]),cimag(s11[i]),delay,freqq2[i]);
    //    {for(i=0;i<np;i++) printf("k %d i %d data %f s11 %f
    //    %f\n",k,i,data[i],creal(s11[i]),cimag(s11[i]));}
    polyfitr(nfit, np, data, mcalc, wtt, dataout, fitf, aarr, bbrr);
    //    for(i=0;i<np;i++) printf("k %d i %d data %f dataout %f diff
    //    %f\n",k,i,data[i],dataout[i],data[i]-dataout[i]);
    for (i = 0; i < npout; i++) {
      sum = 0;
      a = 0;
      for (j = 0; j < nfit; j++) {
        if ((nfit > 16) | (model == 0)) {  //  
          if (j == 0) a = 1;
          if (j % 2) a = cos((freqq1[i] - f0) * 2.0 * PI * (j / 2 + 1) / dfreq);
          if (j % 2 == 0 && j > 0) a = sin((freqq1[i] - f0) * 2.0 * PI * (j / 2) / dfreq);
        } else if ((model == 1) | (model < 0)) {                    // 
          a = pow(log10(freqq1[i] / fcen), j);  //
        } else {
          printf("model must be -1, 0 or 1\n");
          return;
        }

        sum += bbrr[j] * a;
      }
      if (k == 0)
        Taa[i] = sum;
      else
        Taa[i] += sum * I;
    }
    for (i = 0; i < np; i++) {
      sum = 0;
      a = 0;
      for (j = 0; j < nfit; j++) sum += bbrr[j] * fitf[i + j * np];
      if (k == 0)
        Taaa[i] = sum;
      else
        Taaa[i] += sum * I;
    }
  }
  for (i = 0; i < np; i++) Taaa[i] = Taaa[i] * cexp(-2.0 * PI * freqq2[i] * 1e6 * delay * I);
  for (i = 0; i < npout; i++) Taa[i] = Taa[i] * cexp(-2.0 * PI * freqq1[i] * 1e6 * delay * I);
  sum = 0;
  rmsdb = rmsphase = 0;
  for (i = 0; i < np; i++) {
    dbdiff[i] = a = 20.0 * log10(cabs(s11[i]) + 1e-6) - 20.0 * log10(cabs(Taaa[i]) + 1e-6);
    rmsdb += wtt[i] * a * a;
    a = (atan2(cimag(s11[i]), creal(s11[i])) - atan2(cimag(Taaa[i]), creal(Taaa[i]))) * 180.0 / PI;
    if (a > 180.0) a += -360.0;
    if (a < -180.0) a += 360.0;
    pdiff[i] = a;
    rmsphase += wtt[i] * a * a;
    sum += wtt[i];
    if (wtt[i] == 0) pdiff[i] = dbdiff[i] = 0;
  }
  rmsdb = sqrt(rmsdb / sum);
  rmsphase = sqrt(rmsphase / sum);
  if (mode == 0) printf("ant  S11 rms %f dB rms %f deg nfit4 %d %5.1f ns np %d npout %d\n", rmsdb, rmsphase, nfit, delay * 1e9, np, npout);
  if (mode == 1) printf("LNA  S11 rms %f dB rms %f deg nfit3 %d %5.1f ns np %d npout %d\n", rmsdb, rmsphase, nfit, delay * 1e9, np, npout);
  if (mode == 2) printf("hot  S11 rms %f dB rms %f deg nfit2 %d %5.1f ns np %d npout %d\n", rmsdb, rmsphase, nfit, delay * 1e9, np, npout);
  if (mode == 3) printf("amb  S11 rms %f dB rms %f deg nfit2 %d %5.1f ns np %d npout %d\n", rmsdb, rmsphase, nfit, delay * 1e9, np, npout);
  if (mode == 4) printf("cab1 S11 rms %f dB rms %f deg nfit2 %d %5.1f ns np %d npout %d\n", rmsdb, rmsphase, nfit, delay * 1e9, np, npout);
  if (mode == 5) printf("cab2 S11 rms %f dB rms %f deg nfit2 %d %5.1f ns np %d npout %d\n", rmsdb, rmsphase, nfit, delay * 1e9, np, npout);
  if (mode == 6) printf("tmp rms %f dB rms %f deg nfit1 %d np %d npout %d\n", rmsdb, rmsphase, nfit, np, npout);
  if (mode == 7) printf("tmp rms %f dB rms %f deg nfit2 %d np %d npout %d\n", rmsdb, rmsphase, nfit, np, npout);
  if (mode == 8) printf("tmp rms %f dB rms %f deg nfit2 %d np %d npout %d\n", rmsdb, rmsphase, nfit, np, npout);
  if (mode == 9) printf("tmp rms %f dB rms %f deg nfit2 %d np %d npout %d\n", rmsdb, rmsphase, nfit, np, npout);
  if (mode <= 5) plotvna(np, freqq2, s11, wtt, npout, freqq1, Taa, dbdiff, pdiff, mode, rmsdb, rmsphase, nfit, fname);
}

double wavemodel(complex double Ta, complex double Tl, double tsky, double tlnau, double tlnacc, double tlnacs, double tlna) {
  complex double F;
  double t, ph;
  F = sqrt(1.0 - Tl * conj(Tl)) / (1.0 - Ta * Tl);  // F = 1.0/(1.0 - Ta*Tl) gets same results after cal
  ph = atan2(cimag(Ta * F), creal(Ta * F));
  //  ph = atan2(cimag(Ta*F*Tl),creal(Ta*F*Tl));
  t = (tsky * (1.0 - Ta * conj(Ta)) + tlnau * Ta * conj(Ta)) * F * conj(F) + cabs(Ta) * cabs(F) * (tlnacc * cos(ph) + tlnacs * sin(ph)) + tlna;
  return t;
}

double w3p(complex double Ta, complex double Tl, double tsky, double tload, double tamb, double tlnau, double tlnacc, double tlnacs) {
  double pant, pcal, pload;
  pant = wavemodel(Ta, Tl, tsky, tlnau, tlnacc, tlnacs, 0.0);
  pload = wavemodel(0, Tl, tload, tlnau, tlnacc, tlnacs, 0.0);
  pcal = wavemodel(0, Tl, 100.0 + tload, tlnau, tlnacc, tlnacs, 0.0);
  return 100.0 * (pant - pload) / (pcal - pload) + tamb;  // assumes first stage processing uses same tamb and tcal for
                                                          // all spectra
}

double w3pinv(complex double Ta, complex double Tl, double tin, double tload, double tamb, double tlnau, double tlnacc, double tlnacs) {
  complex double F;
  double tout;
  F = 1.0 / (1.0 - Ta * Tl);
  tout = tin - w3p(Ta, Tl, 0.0, tload, tamb, tlnau, tlnacc, tlnacs);
  tout = tout / ((1 - Ta * conj(Ta)) * F * conj(F));
  return tout;
}

double antloss(double freq, double mcalc[], double fitf[], long double aarr[], double bbrr[], double aloss, int mode) {
  double sum, sum2;
  int i, j, nfit, neff;
  double data[100], wt[100];
  double b[10];
  // eff from FEKO 100 to 195 MHz
  //    double eff[] =
  //    {99.9066,99.9254,99.9388,99.9479,99.9534,99.9559,99.9560,99.9538,99.9496,99.9437,99.9361,99.9269,99.9162,99.9040,99.8904,99.8755,99.8592,99.8417,99.8229,99.8032};
  //  eff from FEKO 100 to 200 MHz Blade  90 to 95 99.8965 and 99.9185
  double eff[] = {99.9364, 99.9485, 99.9576, 99.9644, 99.9694, 99.9729, 99.9751, 99.9764, 99.9769, 99.9766, 99.9758,
                  99.9745, 99.9727, 99.9706, 99.9681, 99.9654, 99.9625, 99.9594, 99.9561, 99.9528, 99.9493};
  //    double los[] = {3.708375e+01,  -1.174062e+02,  1.632379e+02,
  //    -1.067407e+02,  2.703088e+01}; double los[] = {4.412460e+01,
  //    -1.632381e+02,  2.202354e+02,  -1.245132e+02,  2.495915e+01};
  double losh[] = {3.470102e+01, -1.335697e+02, 1.824230e+02, -1.000643e+02, 1.874884e+01};  // highband GFhigh
  //    double los[] = {1.564797e+02,  -5.130209e+02,  6.599304e+02,
  //    -3.749992e+02,  7.931530e+01}; // temploss9 - 2.5x2.5 0.1 genblade9e GF
  //    13 2e-3 double los[] = {7.580215e+01,  -2.684529e+02,  3.635722e+02,
  //    -2.154928e+02,  4.726889e+01}; // GF 3.5 1e-3
  double los[] = {5.744029e+01, -1.726088e+02, 2.142321e+02, -1.225013e+02, 2.702699e+01};  // GF 3.5 180deg
  //    double losk[] = {6.410025e+00,  -1.549803e+01,  1.544351e+01,
  //    -7.041051e+00,  1.300656e+00}; double losk[] = {6.321415e+01,
  //    -1.886167e+02,  2.293411e+02,  -1.206021e+02,  2.315841e+01}; // for
  //    low3 memo 290 double losk[] = {1.041479e+01,
  //    -6.031950e+00,  2.307908e+00}; // for low3 3-terms double losk[] =
  //    {-1.179120134559797e+01,  8.816080087652101e+01, -2.263158238997412e+02,
  //    2.823787777684456e+02,  -1.647342754171209e+02,  3.580581311748211e+01};
  double losk[] = {-1.345035060985746e+00, 2.240940678524730e+00, 7.902278433276813e+00, -5.318319232436456e+00};  // loss_low3a_2e-2_rock.txt -nfit 4
                                                                                                                   // -fmin 52 -fmax 110
  double losg5[] = {4.370180753656857e+02, -7.610893859849077e+02, 7.252089244174796e+02, -3.487898288054475e+02,
                    6.949783330081493e+01};  // ground loss
  double losg7[] = {4.846670825246866e+02, -1.023102678080684e+03, 1.303906332843779e+03, -1.000317701465261e+03,
                    4.580417451841646e+02, -1.129872060290666e+02, 1.163765940857081e+01};  // ground loss
  double losg8[] = {3.926509849561377e+00, -1.838394495633041e+01, 3.929777791337121e+01, -4.597706096084568e+01,
                    3.037985609919453e+01, -1.062848140609678e+01, 1.531802534710291e+00};  // mesh loss

  if (mode == 1) {  // only use antloss for lmode = 1
    nfit = 7;
    neff = 20 + 1;
    for (j = 0; j < nfit; j++) {
      for (i = 0; i < neff; i++) {
        data[i] = eff[i];
        wt[i] = 1;
        fitf[i + j * neff] = pow(100 + i * 5 - 150.0, j);
      }
    }
    polyfitr(nfit, neff, data, mcalc, wt, data, fitf, aarr, bbrr);
    for (j = 0; j < nfit; j++) b[j] = bbrr[j];
    sum = 0;
    for (j = 0; j < nfit; j++) sum += b[j] * pow(freq - 150.0, j);
    sum = sum * 1e-2;
    nfit = 5;
    sum2 = 0;
    for (j = 0; j < nfit; j++) sum2 += losh[j] * pow(freq / 150.0, j);  // when using GFhigh
    sum2 = 0.5e-2 * 300.0;                                              // assume 0.5% groundloss
    return sum * (1.0 - sum2 / 300.0);                                  // antloss plus ground loss
  }
  if (mode == 5) {  // for lmode = 5 original ground plane
    nfit = 5;
    sum = 0;
    for (j = 0; j < nfit; j++) sum += los[j] * pow(freq * 2.0 / 150.0, j);  // when using GF
    //     for(j=0;j<nfit;j++) sum += (los[j]/4.0)*pow(freq * 2.0 / 150.0,j); //
    //     divide by 4 to get 5x5m printf("freq %f loss %f\n",freq,sum);
    return (1.0 - sum / 300.0);
  }
  if (mode == 7 || mode == 9) {  // 0.5 percent
    return (1.0 - 0.5e-2);
  }
  if (mode == 8) {
    sum = 0;
    nfit = 4;  // low3
    for (j = 0; j < nfit; j++) sum += losk[j] * pow(freq / 75.0, j);
    //     printf("freq %f loss %f\n",freq,sum/300.0);
    return (1.0 - sum / 300.0);
  }
  if (mode == 20) {
    sum = 0;
    nfit = 5;  // EDGES-3 ground loss
    for (j = 0; j < nfit; j++) sum += losg5[j] * pow(freq / 75.0, j);
    //     printf("freq %f loss %f\n",freq,sum/300.0);
    return (1.0 - sum / 300.0);
  }
  if (mode == 21) {
    sum = 0;
    nfit = 7;  // EDGES-3 ground loss
    for (j = 0; j < nfit; j++) sum += losg7[j] * pow(freq / 75.0, j);
    return (1.0 - sum / 300.0);
  }
  if (mode == 22) {
    sum = 0;
    nfit = 7;  // EDGES-3 mesh loss
    for (j = 0; j < nfit; j++) sum += losg8[j] * pow(freq / 75.0, j);
    return (1.0 - sum / 300.0);
  }
  if (mode == 100) {
    return (pow(10.0,
                -0.1 * aloss * sqrt(freq / 100.0)));  // loss in dB independent of S11
  }
  return 1;
}

double loss(complex double s11a, double atten)  // simple attenuation
{
  complex double T, ss11, ss12, ss21, ss22;
  double L, La;
  La = pow(10.0, -0.1 * atten);  // atten in dB
  ss11 = ss22 = 0;
  ss21 = sqrt(La);
  ss12 = sqrt(La);
  T = (s11a - ss11) / (ss12 * ss21 - ss11 * ss22 + ss22 * s11a);  // from memo 132
  L = ss21 * conj(ss21) * (1 - T * conj(T)) / ((1 - s11a * conj(s11a)) * (1 - ss22 * T) * (1 - conj(ss22 * T)));
  //   if(L < 0.1) printf("error check loss L=%f\n",L);
  return L;
}

double lossmodel(double freq, complex double s11a, int mode) {
  complex double s11, s22, s12, s21, ta11, ta12, ta21, ta22, tb11, tb12, tb21, tb22, t11, t12, t21, t22, T;
  double del, del4, L2, ttest, ttest2;
  double s11r[] = {0.000105, 0.006557, -0.001168, -0.000260};
  double s11i[] = {0.000907, 0.003658, -0.005924, 0.001011};
  double s12r[] = {0.995152, 0.011466, -0.326861, 0.043317};
  double s12i[] = {0.000360, -0.767049, 0.034407, 0.053560};
  double s22r[] = {0.004139, -0.011762, 0.021795, -0.009117};
  double s22i[] = {0.000500, 0.001785, 0.000988, -0.002854};

  //  double ss11r[] = {0.000367, 0.005528, -0.000934, -0.000061};
  //  double ss11i[] = {0.000548, 0.004290, -0.005296, 0.001079};
  //  double ss12r[] = {0.998282, 0.002708, -0.317147, 0.041610};
  //  double ss12i[] = {-0.000086, -0.764256, 0.038001, 0.050273};
  //  double ss22r[] = {0.000221, 0.005526, -0.002201, -0.000244};
  //  double ss22i[] = {0.000466, 0.002232, -0.006258, 0.001766};

  double ss11r[] = {0.000622, 0.004519, 0.000199, -0.000442};  // lossfit2m2a 40 - 200 using
                                                               // semi_rigid_s_parameters_WITH_HEADER.txt
  double ss11i[] = {0.000856, 0.003064, -0.003920, 0.000617};
  double ss12r[] = {0.994266, 0.017943, -0.333811, 0.047113};
  double ss12i[] = {0.002412, -0.773510, 0.047980, 0.047009};
  double ss22r[] = {0.000627, 0.004000, -0.000539, -0.000792};
  double ss22i[] = {0.000832, 0.000768, -0.004606, 0.001209};

  int m;
  // mode 0 = hot load  5 = low band 5 or 6 or 7 = lowband
  if (mode == 3) return lossmodel2(freq, s11a);
  if (mode >= 20) return 1;  // EDGES-3
  ttest = ttest2 = 1;        // ttest = 0.5;
  if (mode == 9) ttest2 = 0.8;
  del4 = 1e-10;                                                                                    // lowband Fairview SC3792
  if (mode == 1) del4 = 1e-10;                                                                     // needs to be checked
  if (mode == 11) del4 = 0.9 * 2.54e-2 / 3e08;                                                     // midband
  del = 21.90 * 2.54e-2 / 3e08;                                                                    // balun  high band 21"
  if (mode == 5 || mode == 6 || mode == 7 || mode == 8 || mode == 9) del = 43.6 * 2.54e-2 / 3e08;  // balun lowband 43"

  if (mode == 10) del = 35 * 2.54e-2 / 3e08;    // balun midband 35"
  if (mode == 11) del = 33.5 * 2.54e-2 / 3e08;  // balun midband 33.5"
  //   if(mode==4) {del4=2.5*2.54e-2/3e8; del=0.51/3e08; ttest=0.2;} // balun +
  //   bends for active coldload test if(mode==0) {del=0; del4=3.03e-10;}  //
  //   for hotload  - fit to semi from raul
  if (mode == 5 || mode == 6 || mode == 7 || mode == 8 || mode == 9 || mode == 10 || mode == 11) {
    cabl2(freq, del4, 0, 8, &s11, &s12, ttest);
    s21 = s12;
    s22 = s11;
  }  // SC3792
  else {
    cabl2(freq, del4, 0, 9, &s11, &s12, ttest);
    s21 = s12;
    s22 = s11;
  }                 // bend or stainless 0.141
  if (mode == 2) {  // fit to /home/aeer/px14/mro/data2/semi_rigid_S_parameters.txt - as
                    // example
    s11 = 0;
    for (m = 0; m < 4; m++) s11 += s11r[m] * pow(freq / 150.0, m) + s11i[m] * pow(freq / 150.0, m) * I;
    s12 = 0;
    for (m = 0; m < 4; m++) s12 += s12r[m] * pow(freq / 150.0, m) + s12i[m] * pow(freq / 150.0, m) * I;
    s22 = 0;
    for (m = 0; m < 4; m++) s22 += s22r[m] * pow(freq / 150.0, m) + s22i[m] * pow(freq / 150.0, m) * I;
    s21 = s12 = csqrt(s12);
    T = (s11a - s11) / (s12 * s21 - s11 * s22 + s22 * s11a);  // from memo 132 - need to check
                                                              // Raul's convention s11 vs s22
    return cabs(s12 * s21) * (1 - T * conj(T)) / ((1 - s11a * conj(s11a)) * (1 - s22 * T) * (1 - conj(s22 * T)));
  }
  if (mode == 0) {  // fit to /home/aeer/px14/mro/data4/semi_rigid_s_parameters.txt - as
                    // example lossfit2m2.c
    s11 = 0;
    for (m = 0; m < 4; m++) s11 += ss11r[m] * pow(freq / 150.0, m) + ss11i[m] * pow(freq / 150.0, m) * I;
    s12 = 0;
    for (m = 0; m < 4; m++) s12 += ss12r[m] * pow(freq / 150.0, m) + ss12i[m] * pow(freq / 150.0, m) * I;
    s22 = 0;
    for (m = 0; m < 4; m++) s22 += ss22r[m] * pow(freq / 150.0, m) + ss22i[m] * pow(freq / 150.0, m) * I;
    s21 = s12 = csqrt(s12);
    T = (s11a - s11) / (s12 * s21 - s11 * s22 + s22 * s11a);  // from memo 132 - need to check
                                                              // Raul's convention s11 vs s22
    return cabs(s12 * s21) * (1 - T * conj(T)) / ((1 - s11a * conj(s11a)) * (1 - s22 * T) * (1 - conj(s22 * T)));
  }
  if (mode == 12) {  // loss for Molex cable
    del4 = 10 * 12 * 2.54e-2 / 3e08;
    cabl2(freq, del4, 0, 12, &s11, &s12, ttest);
    s21 = s12;
    s22 = s11;
    T = (s11a - s11) / (s12 * s21 - s11 * s22 + s22 * s11a);  // from memo 132 - need to check
                                                              // Raul's convention s11 vs s22
    return cabs(s12 * s21) * (1 - T * conj(T)) / ((1 - s11a * conj(s11a)) * (1 - s22 * T) * (1 - conj(s22 * T)));
  }
  ta11 = -(s11 * s22 - s12 * s21) / s21;
  ta12 = s11 / s21;
  ta21 = -s22 / s21;
  ta22 = 1 / s21;  // close to ref. plane
  if (mode == 1) {
    cabl2(freq, del, 0, 0, &s11, &s12, 1.0);
    s21 = s12;
    s22 = s11;
  }  // high band balun tube
  if (mode == 5 || mode == 6 || mode == 7 || mode == 8 || mode == 9) {
    cabl2(freq, del, 0, 5, &s11, &s12, ttest2);
    s21 = s12;
    s22 = s11;
  }  // low band balun tube
  if (mode == 10 || mode == 11) {
    cabl2(freq, del, 0, 21, &s11, &s12, ttest2);
    s21 = s12;
    s22 = s11;
  }  // mid band balun tube
  tb11 = -(s11 * s22 - s12 * s21) / s21;
  tb12 = s11 / s21;
  tb21 = -s22 / s21;
  tb22 = 1 / s21;  // T matrix from s-parms closer to antenna
  t11 = ta11 * tb11 + ta12 * tb21;
  t12 = ta11 * tb12 + ta12 * tb22;
  t21 = ta21 * tb11 + ta22 * tb21;
  t22 = ta21 * tb12 + ta22 * tb22;  // A*B
  //   t11=tb11*ta11+tb12*ta21; t12=tb11*ta12+tb12*ta22;
  //   t21=tb21*ta11+tb22*ta21; t22=tb21*ta12+tb22*ta22; //B*A
  s11 = t12 / t22;
  s12 = (t11 * t22 - t12 * t21) / t22;
  s21 = 1 / t22;
  s22 = -t21 / t22;
  T = (s11a - s11) / (s12 * s21 - s11 * s22 + s22 * s11a);  // from memo 132
  L2 = cabs(s12 * s21) * (1 - T * conj(T)) / ((1 - s11a * conj(s11a)) * (1 - s22 * T) * (1 - conj(s22 * T)));
  //   T=cabl2(freq,del,T,5,&s11,&s12,1.0); T =
  //   cabl2(freq,del4,T,8,&s11,&s12,ttest); printf("check s11a %f %f  %f
  //   %f\n",creal(s11a),cimag(s11a),creal(T),cimag(T));
  return L2;
}

double rigloss(complex double s11a, complex double s11, complex double s12, complex double s22) {
  complex double T;
  T = (s11a - s11) / (s12 - s11 * s22 + s22 * s11a);  // from memo 132  s12 = s12*s21
  return cabs(s12) * (1 - T * conj(T)) / ((1 - s11a * conj(s11a)) * (1 - s22 * T) * (1 - conj(s22 * T)));
}

double lossmodel2(double freq, complex double s11a) {
  int m;
  double L;
  complex double T, s11, s12, s21, s22;
  // from fit to s11 data using open,short,load

  double s11r[] = {-0.021797, 0.134691, -0.162748, 0.053239};
  double s11i[] = {0.004280, 0.025673, -0.089287, 0.051092};
  double s12r[] = {0.993659, 0.012526, -0.376041, 0.053958};
  double s12i[] = {-0.000998, -0.818758, 0.039535, 0.065117};
  double s22r[] = {-0.015862, 0.110943, -0.144476, 0.053251};
  double s22i[] = {0.012254, -0.018995, -0.016018, 0.020065};

  //    double s11r[] = {-0.020739,  0.128536, -0.154593, 0.050743};
  //    double s11i[] = {0.004310,  0.024908, -0.084116, 0.048021};
  //    double s12r[] = {0.992806,  0.013482, -0.376738, 0.054316};
  //    double s12i[] = {-0.001003,  -0.818058, 0.039819, 0.064806};
  //    double s22r[] = {-0.014396,  0.099783, -0.126599, 0.045170};
  //    double s22i[] = {0.007887,  -0.001775, -0.034834, 0.025757};

  s11 = 0;
  for (m = 0; m < 4; m++) s11 += s11r[m] * pow(freq / 150.0, m) + s11i[m] * pow(freq / 150.0, m) * I;
  s12 = 0;
  for (m = 0; m < 4; m++) s12 += s12r[m] * pow(freq / 150.0, m) + s12i[m] * pow(freq / 150.0, m) * I;
  s22 = 0;
  for (m = 0; m < 4; m++) s22 += s22r[m] * pow(freq / 150.0, m) + s22i[m] * pow(freq / 150.0, m) * I;
  s21 = s12 = csqrt(s12 * cexp(-2 * PI * freq * 1e6 * 0.5 / 3e8 * 2 * I));
  T = (s11a - s11) / (s12 * s21 - s11 * s22 + s22 * s11a);
  L = s21 * conj(s21) * (1 - T * conj(T)) / ((1 - s11a * conj(s11a)) * (1 - s22 * T) * (1 - conj(s22 * T)));
  //     printf("freq %f s11 %f %f %f s12 %f %f %f s22 %f %f %f Loss
  //     %f\n",freqs[i],creal(s11),cimag(s11),cabs(s11),creal(s12*s21),cimag(s12*s21),cabs(s12*s21),creal(s22),cimag(s22),cabs(s22),L);
  return L;
}

complex double cabl2(double freq, double delay, complex double Tin, int mode, complex double *ss11, complex double *ss12, double ttest) {
  complex double T, Zcab, g, s11, s12, s21, s22, Vin, Iin, Vout, VVin, Z;
  double a, b, d, d2, diel, R, C, L, La, Lb, disp, G;
  // ttest = 0.5;  // best fit from day 2015_023
  {
    b = 0.37 * 2.54e-2 * 0.5;
    a = (5.0 / 32.0) * 2.54e-2 * 0.5;
    diel = 1.07;                                               // balun tube
    d2 = sqrt(1.0 / (PI * 4.0 * PI * 1e-7 * 5.96e07 * 0.29));  // skin depth at 1 Hz for brass
    d = sqrt(1.0 / (PI * 4.0 * PI * 1e-7 * 5.96e07));          // skin depth at 1 Hz for copper
  }
  if (mode == 5) {
    b = 0.75 * 2.54e-2 * 0.5;
    a = (5.0 / 16.0) * 2.54e-2 * 0.5;
    diel = 1.07;                                               // balun tube lowband
    d2 = sqrt(1.0 / (PI * 4.0 * PI * 1e-7 * 5.96e07 * 0.29));  // skin depth at 1 Hz for brass
    d = sqrt(1.0 / (PI * 4.0 * PI * 1e-7 * 5.96e07 * ttest));  // skin depth at 1 Hz for copper
  }
  if (mode == 8) {
    b = 0.161 * 2.54e-2 * 0.5;
    a = 0.05 * 2.54e-2 * 0.5;
    diel = 2.05;                                                        // SC3792 connector - new dimensions from Fairviwe 8 Dec 15
    d2 = sqrt(1.0 / (PI * 4.0 * PI * 1e-7 * 5.96e07 * 0.024 * ttest));  // for Stainless
    d = sqrt(1.0 / (PI * 4.0 * PI * 1e-7 * 5.96e07 * 0.24 * ttest));    // skin depth at 1 Hz for brass - might be less 0.24
  }
  if (mode == 9) {
    b = 0.16 * 2.54e-2 * 0.5;
    a = 0.05 * 2.54e-2 * 0.5;
    diel = 2.05;                                                        // SMA connector
    d2 = sqrt(1.0 / (PI * 4.0 * PI * 1e-7 * 5.96e07 * 0.024 * ttest));  // for Stainless
    d = sqrt(1.0 / (PI * 4.0 * PI * 1e-7 * 5.96e07 * 0.20 * ttest));    // skin depth at 1 Hz note change
  }
  if (mode == 12) {
    b = 0.1515 * 2.54e-2 * 0.5;
    a = 0.0453 * 2.54e-2 * 0.5;
    diel = 2.1;                                         // Molex WM10479-nd FEP molex 100054008 could only find inner =
                                                        // 0.0453 inch twicked outer to get 50 ohms
    d2 = sqrt(1.0 / (PI * 4.0 * PI * 1e-7 * 5.96e07));  // for silver plated copper
    d = sqrt(1.0 / (PI * 4.0 * PI * 1e-7 * 5.96e07));   // for silver plated copper
  }
  if (mode == 21) {
    b = 1.25 * 2.54e-2 * 0.5;
    a = (16.0 / 32.0) * 2.54e-2 * 0.5;
    diel = 1.2;                                                // balun tube midband  8859k332    7782T333
    d2 = sqrt(1.0 / (PI * 4.0 * PI * 1e-7 * 5.96e07 * 0.29));  // skin depth at 1 Hz for brass
    d = sqrt(1.0 / (PI * 4.0 * PI * 1e-7 * 5.96e07 * ttest));  // skin depth at 1 Hz for copper
  }
  L = (4.0 * PI * 1e-7 / (2.0 * PI)) * log(b / a);
  C = 2.0 * PI * 8.854e-12 * diel / log(b / a);

  La = 4.0 * PI * 1e-7 * d / (4.0 * PI * a);
  Lb = 4.0 * PI * 1e-7 * d2 / (4.0 * PI * b);
  disp = (La + Lb) / L;
  R = 2.0 * PI * L * disp * sqrt(freq * 1e6);
  L = L * (1.0 + disp / sqrt(freq * 1e6));
  G = 0;
  if (diel > 1.2) G = 2.0 * PI * C * freq * 1e6 * 2e-4;
  Zcab = csqrt((I * 2 * PI * freq * 1e6 * L + R) / (I * 2 * PI * freq * 1e6 * C + G));
  g = csqrt((I * 2 * PI * freq * 1e6 * L + R) * (I * 2 * PI * freq * 1e6 * C + G));

  T = (50.0 - Zcab) / (50.0 + Zcab);
  Vin = (cexp(+g * delay * 3e08) + T * cexp(-g * delay * 3e08));
  Iin = (cexp(+g * delay * 3e08) - T * cexp(-g * delay * 3e08)) / Zcab;
  Vout = (1 + T);  // Iout = (1 - T)/Zcab;
  s11 = s22 = ((Vin / Iin) - 50) / ((Vin / Iin) + 50);
  VVin = Vin + 50.0 * Iin;
  s12 = s21 = (2 * Vout / VVin);
  *ss11 = s11;
  *ss12 = s12;

  Z = 50.0 * (1 + Tin) / (1 - Tin);
  T = (Z - Zcab) / (Z + Zcab);
  T = T * cexp(-g * 2 * delay * 3e08);
  Z = Zcab * (1 + T) / (1 - T);
  T = (Z - 50.0) / (Z + 50.0);
  return T;
}

/*
Lesurf's approx. maybe slightly better 6 Nov 15
complex double cabl2(double freq, double delay, complex double Tin, int
mode,complex double *ss11, complex double *ss12, double ttest) { complex double
T,Zcab,g,s11,s12,s21,s22,Vin,Iin,Vout,Iout,VVin,Z; double a,b,d,diel,R,C,L,G;
  double q,Rs,R0,Ra,Xa,s,s1,s2;
  // ttest = 0.5;  // best fit from day 2015_023
              {  b=0.37*2.54e-2*0.5; a=(5.0/32.0)*2.54e-2*0.5; diel = 1.07;  //
balun tube s2 = 5.96e07*0.29;  // skin depth at 1 Hz for brass s1 = 5.96e07;  //
skin depth at 1 Hz for copper
              }
  if(mode==5) {  b=0.75*2.54e-2*0.5; a=(5.0/16.0)*2.54e-2*0.5; diel = 1.07;  //
balun tube lowband s2 = 5.96e07*0.29;  // skin depth at 1 Hz for brass s1
= 5.96e07;  // skin depth at 1 Hz for copper
              }
  if(mode==8) { b=0.304*2.54e-2*0.5; a=0.086*2.54e-2*0.5; diel = 2.05;  //
SC3792 connector s2 = 5.96e07*0.024*ttest;   // for Stainless s1
= 5.96e07*0.24*ttest;  // skin depth at 1 Hz for brass - might be less 0.24
              }
  if(mode==9) { b=0.16*2.54e-2*0.5; a=0.05*2.54e-2*0.5; diel = 2.05;  // SMA
connector s2 = 5.96e07*0.024*ttest;   // for Stainless s1 = 5.96e07*0.05*ttest;
// skin depth at 1 Hz for brass - might be less 0.24
              }
  L = (4.0*PI*1e-7/(2.0*PI))*log(b/a);
  C=2.0*PI*8.854e-12*diel/log(b/a);

  s = s1;
  d = sqrt(1.0/(PI*4.0*PI*1e-7*s));
  q = sqrt(2.0)*a*sqrt(freq*1e6)/d;
  Rs = sqrt(PI*freq*1e6*4.0*PI*1e-7/s);
  R0 = 2.0*Rs/(sqrt(2.0)*PI*a*q);
  Ra = R0*(q/(sqrt(2.0)*2) + 0.26);   // Lesurf's approx
  Xa = R0*(q/(sqrt(2.0)*2) - 0.02);

  s = s2;
  d = sqrt(1.0/(PI*4.0*PI*1e-7*s));
  q = sqrt(2.0)*b*sqrt(freq*1e6)/d;
  Rs = sqrt(PI*freq*1e6*4.0*PI*1e-7/s);
  R0 = 2.0*Rs/(sqrt(2.0)*PI*b*q);
  Ra += R0*(q/(sqrt(2.0)*2) + 0.26);
  Xa += R0*(q/(sqrt(2.0)*2) - 0.02);

  R = Ra;
  L += Xa/(2*PI*freq*1e6);

  G=0;
  if(diel > 1.2) G=2.0*PI*C*freq*1e6*2e-4;
  Zcab = csqrt((I*2*PI*freq*1e6*L+R)/(I*2*PI*freq*1e6*C+G));
  g = csqrt((I*2*PI*freq*1e6*L+R)*(I*2*PI*freq*1e6*C+G));

  T = (50.0-Zcab)/(50.0+Zcab);
  Vin = (cexp(+g*delay*3e08) + T*cexp(-g*delay*3e08));
  Iin = (cexp(+g*delay*3e08) - T*cexp(-g*delay*3e08))/Zcab;
  Vout = (1 + T); Iout = (1 - T)/Zcab;
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
*/

double lossinv(double tant, double tamb, double L) { return (tant - tamb * (1 - L)) / L; }

double beamcorr(double freq, int bfit, double ang, double secs, double mcalc[], double fitf[], long double aarr[], double bbrr[], int skymode,
                int integ, double *lst, double bb[], double *bspac, int cmb, double *freqref, double fbstart, double fbstop, int site, int map) {
  char name[255], buf[32768], *p;
  FILE *file3;
  static double azel[360 * 91 * NBEAM], data[NBEAM], data2[NBEAM], wt[NBEAM], fsum[NBEAM], fsum1[NBEAM], fsum2[NBEAM], poww[NBEAM], mapgal[512][1024];
  static double rraa[512][1024], sindec[512][1024], cosdec[512][1024];
  static double map45[512][1024], spind[512][1024];
  double k, sum, amp, amp150, raa, dec, lat, lon, gstt, sunra, sundec, sunaz, sunel, glat, glon, wsum, azz, el, sang, mp, dp, gp, opac, telec, rms,
      wb, aind, wsum2;
  double ssecs, fr, cmbb, av45, av408, nav;
  int i, j, n, nn, az, frq, frq150, iaz, iel, m, frqst, frqspac, nbeam, yr, dy, hr, mn, sc, ii, jj, iii, jjj;
  double cost, sint, zb, antght, d, max;
  double sinlat, coslat;
  frqst = 0;
  frqspac = 0;
  n = 0;
  max = 0;
  if (integ > 0) {
    if (cmb & 2)
      cmbb = 2.725;
    else
      cmbb = 0;
    i = -1;
    sprintf(name, "azelq.txt");
    if ((file3 = fopen(name, "r")) == NULL) {
      printf("%s error\n", name);
      return 0;
    }
    while (fgets(buf, 32768, file3) != 0) {
      if (strstr(buf, "#FREQ")) {
        sscanf(buf, "%*s %d", &frq);
        i++;
        if (i == 0 && fbstart < frq) frqst = frq;
        if (frq == fbstart) {
          frqst = frq;
          i = 0;
        }
        if (frq == fbstop) n = i;
        if (i == 1) frqspac = frq - frqst;
      }
      if (buf[0] == 'a' && i >= 0 && i < NBEAM) {
        p = &buf[3];
        j = 0;
        az = 0;
        az = strtod(p, &p) + ang;
        if (az < 0) az += 360;
        if (az < 0) az += 360;
        if (az >= 360) az -= 360;

        while (j < 91 && p && *p != '\n') {
          azel[(az + j * 360) * NBEAM + i] = pow(10.0, 0.1 * strtod(p, &p));
          //                           printf("az %d el %d fr %d amp
          //                           %f\n",az,j,frqst+i*frqspac,azel[(az+j*360)*NBEAM+i]);

          if (frqst < 70.0)
            *freqref = 75.0;
          else
            *freqref = 150.0;
          if (skymode >= 512) {
            el = j;
            cost = cos(el * PI / 180.0) * cos((az + ang) * PI / 180.0);
            zb = sin(el * PI / 180.0);
            sint = sqrt(1.0 - cost * cost);
            //  if(fabs(sint)>1e-6 && el>=0.0) amp=(cos(0.5*PI*cost))/sint;
            if (fabs(sint) > 1e-6 && el >= 0.0) {
              amp = PI * ((frqst + i * frqspac) / 150.0) / 2.0;  // half-wavelength at 150 MHz dipole
              amp = (cos(amp * cost) - cos(amp)) / sint;
            } else
              amp = 1e-6;
            antght = 0.25;  // height in wavelengths at 150 MHz
            if (zb > 0.0)
              amp = amp * sin(PI * 2.0 * zb * antght * (frqst + i * frqspac) / 150.0);
            else
              amp = 1e-6;
            // if((az==0 || az==90) && el==60)     printf("az %d el %d fr %d
            // amp1 %f amp2 %f ratio
            // %f\n",az,j,frqst+i*frqspac,azel[(az+j*360)*NBEAM+i],amp*amp,azel[(az+j*360)*NBEAM+i]/(amp*amp));
            azel[(az + j * 360) * NBEAM + i] = amp * amp;
          }
          j++;
        }
      }
    }
    fclose(file3);
    if (i < NBEAM + 1)
      nbeam = i + 1;
    else
      return 0;
    if (n) nbeam = n + 1;
    // for(i=0;i<nbeam;i++) printf("i %d beam %f %f %f
    // %f\n",i,azel[(0+0*360)*NBEAM+i],azel[(1+0*360)*NBEAM+i],azel[(0+1*360)*NBEAM+i],azel[(359+0*360)*NBEAM+i]);
    for (i = 0; i < 512; i++)
      for (j = 0; j < 1024; j++) mapgal[i][j] = -1;
    if (map == 0 || map == 2) {
      if ((file3 = fopen("408-all-noh", "r")) == NULL) { return 0; }
    }  // map 0= use Haslam 1= use Guzman 2=Haslam + Guzman 3=Guzman + Haslam
    if (map == 1 || map == 3) {
      if ((file3 = fopen("/home/aeer/fits/45mhz.txt", "r")) == NULL) { return 0; }
    }  // use Guzmin map
    i = j = 0;
    while (fgets(buf, 32768, file3) != 0) {
      p = buf;
      j = 0;
      k = -1;
      while (j < 1024 && p && sscanf(p, "%lf", &k) == 1) {
        while (*p == ' ') p++;
        if (*p && p) p = strchr(p, ' ');
        if (cmb & 4 && (map == 0 || map == 2)) k += -5;                // offset correction for Haslam 408 map
        if (cmb & 4 && (map == 1 || map == 3)) k = k * 1.076 - 160.0;  // Guzman map correction
        if (cmb & 8) {
          if (i == 256 && j == 512)
            k = 2.3e6;
          else
            k = 18;
        }  // simulate point source
        if (cmb & 64)
          k = exp(-4.0 * log(2.0) * ((i - 256) * (i - 256) + (j - 512) * (j - 512)) / 810.0) * 2.6e3 +
              18;  // 1 unit = 0.3516 deg 28 to simulate source 10 FWHM deg wide
        mapgal[i][j] = k;
        // printf("i %d j %d k %6.0f\n",i,j,k);
        j++;
      }
      i++;
    }
    fclose(file3);
    printf("MAP: %d, SITE %d\n", map, site);
    if (map >= 2) {
      if (map == 2)
        if ((file3 = fopen("/home/aeer/fits/45mhz.txt", "r")) == NULL) { return 0; }
      if (map == 3)
        if ((file3 = fopen("408-all-noh", "r")) == NULL) { return 0; }

      i = j = 0;

      while (fgets(buf, 32768, file3) != 0) {
        p = buf;
        j = 0;
        k = -1;
        while (j < 1024 && p && sscanf(p, "%lf", &k) == 1) {
          while (*p == ' ') p++;
          if (*p && p) p = strchr(p, ' ');
          map45[i][j] = k;
          j++;
        }
        i++;
      }
      fclose(file3);
      for (i = 0; i < 512; i++) {
        for (j = 0; j < 1024; j++) {
          av45 = av408 = nav = 0;
          for (ii = -3; ii <= 3; ii++) {
            for (jj = -3; jj <= 3; jj++) {
              iii = i + ii;
              jjj = j + jj;
              if (iii >= 0 && iii < 512 && jjj >= 0 && jjj < 1024) {
                if (map == 2) av45 += map45[iii][jjj];
                if (map == 2) av408 += mapgal[iii][jjj];
                if (map == 3) av45 += mapgal[iii][jjj];
                if (map == 3) av408 += map45[iii][jjj];
                nav++;
              }
            }
          }
          spind[i][j] = log((av45 / nav - cmbb) / (av408 / nav - cmbb)) / log(408.0 / 45.0);
        }
      }
    }
    telec = 1e3;
    lon = 116.5 * PI / 180.0;  // Boolardy
    lat = -26.7 * PI / 180.0;  // EDGES -26.72 116.61
    if (site == 11) {
      lon = 117.5 * PI / 180.0;
      lat = -27.7 * PI / 180.0;
    }
    if (site == 12) {
      lon = 117.5 * PI / 180.0;
      lat = -26.7 * PI / 180.0;
    }
    if (site == 13) {
      lon = 117.5 * PI / 180.0;
      lat = -25.7 * PI / 180.0;
    }
    if (site == 14) {
      lon = 116.5 * PI / 180.0;
      lat = -27.7 * PI / 180.0;
    }
    if (site == 15) {
      lon = 116.5 * PI / 180.0;
      lat = -26.7 * PI / 180.0;
    }
    if (site == 16) {
      lon = 116.5 * PI / 180.0;
      lat = -25.7 * PI / 180.0;
    }
    if (site == 17) {
      lon = 115.5 * PI / 180.0;
      lat = -27.7 * PI / 180.0;
    }
    if (site == 18) {
      lon = 115.5 * PI / 180.0;
      lat = -26.7 * PI / 180.0;
    }
    if (site == 19) {
      lon = 115.5 * PI / 180.0;
      lat = -25.7 * PI / 180.0;
    }
    if (site == 1) {
      lon = 116.5 * PI / 180.0;
      lat = 42.417 * PI / 180.0;
    }  // Oregon for simulations
    if (site == 2) {
      lon = -119.05 * PI / 180.0;
      lat = 42.417 * PI / 180.0;
    }  // Oregon
    if (site == 3) {
      lon = -90.76 * PI / 180.0;
      lat = 79.433 * PI / 180.0;
    }  // Baffin Island
    if (site == 6) {
      lon = 116.5 * PI / 180.0;
      lat = 79.433 * PI / 180.0;
    }  // Baffin Island for simulation
    if (site == 4) {
      lon = 21.9889 * PI / 180.0;
      lat = -30.969 * PI / 180.0;
    }  // South Africa
    if (site == 5) {
      lon = 116.5 * PI / 180.0;
      lat = 13.0 * PI / 180.0;
    }  // India for simulation
    if (site == 100) {
      lon = 118.5 * PI / 180.0;
      lat = -24.7 * PI / 180.0;
    }  // for test offset from MRO
    if (site > 100) {
      lon += ((site % 10) - 5) * PI / 180.0;
      lat += ((site / 10) - 15) * PI / 180.0;
    }  // 155 = no offset 165 = +1deg lat 156 = +1deg long
    sinlat = sin(lat);
    coslat = cos(lat);
    for (m = 0; m < nbeam; m++) fsum[m] = fsum1[m] = fsum2[m] = 0;
    frq150 = 0;
    k = 1e6;
    for (frq = 0; frq < nbeam; frq++)
      if (fabs(frqst + frq * frqspac - (*freqref)) < k) {
        k = fabs(frqst + frq * frqspac - (*freqref));
        frq150 = frq;
      }  // find closest
    printf("REFERENCE FREQ: %f %d\n", frqst + frq150 * frqspac, frq150);
    for (i = 0; i < 512; i++) {
      glat = (i - 256.0) * 90.0 / 256.0;
      for (j = 0; j < 1024; j++) {
        glon = -(j - 512.0) * 180.0 / 512.0;
        GalactictoRadec(glat, glon, &raa, &dec);
        rraa[i][j] = raa;
        sindec[i][j] = sin(dec);
        cosdec[i][j] = cos(dec);
      }
    }
    opac = 1.8 / 500.0;
    nn = integ / 1800 + 1;  // averaging over 3600 sec makes little difference to result
                            //    nn = integ/900 + 1; // made a small difference
    d = (double)integ / (double)nn;
    printf("USING secs = %f, dsec = %f, nn=%d\n", secs, d, nn);

    for (n = 0; n < nn; n++) {
      ssecs = secs + d / 2.0 + n * d - ((double)integ) / 2.0;
      sunradec(ssecs, &sunra, &sundec);
      gstt = gst(ssecs);
      radec_azel(gstt - sunra + lon, sundec, lat, &sunaz, &sunel);
      sunaz = sunaz * 180.0 / PI;
      sunel = sunel * 180.0 / PI;
      *lst = (gstt + lon) * 12.0 / PI;
      if (*lst < 0) *lst += 24.0;
      if (*lst > 24.0) *lst -= 24.0;
      toyrday(ssecs, &yr, &dy, &hr, &mn, &sc);
      printf(
          "Spsecs %.3f %4d:%03d:%02d:%02d:%02d sunaz %6.2f sunel %6.2f lst "
          "%.6f gha %6.2f\n",
          ssecs, yr, dy, hr, mn, sc, sunaz, sunel, *lst, *lst - (17.0 + 45.67 / 60.0));
      
      sprintf(name, "skymap%d.txt",n);
      file3 = fopen(name, "w");
      fprintf(file3, "az el wsum amp amp150 wsum2 iaz iel sang\n");


      
      for (i = 0; i < 512; i += 1) {
        glat = (i - 256.0) * 90.0 / 256.0;
        sang = (180.0 / 512.0) * cos(glat * PI / 180.0) * (360.0 / 1024.0) * PI * PI / (180.0 * 180.0);
        for (frq = 0; frq < nbeam; frq++) {
          fr = frqst + frq * frqspac;
          if (skymode & 1)
            dp = 0.05;
          else
            dp = 0;                         // spectral index variation
          if (fabs(glat) > 20.0) dp = -dp;  // sign changed nov18
          if (skymode & 2)
            gp = -0.12 * log(fr / 150.0);  // derivative Angelica's gamma
          else
            gp = 0.0;
          
          if ((map == 0) && skymode & 128) poww[frq] = pow(fr / 408.0, -2.5 + dp + gp);
          if ((map == 1) && skymode & 128) poww[frq] = pow(fr / 45.0, -2.5 + dp + gp);
        }
        for (j = 0; j < 1024; j += 1) {
          radec_azel2(gstt - rraa[i][j] + lon, sindec[i][j], cosdec[i][j], sinlat, coslat, &azz, &el);
          azz = azz * 180.0 / PI;
          el = el * 180.0 / PI;
          //  phi=el*PI/180.0; R=6356.0; h=300.0;  opac =
          //  (1.8/500.0)*(R+h)/sqrt(R*R*sin(phi)*sin(phi)+2*R*h+h*h); // from
          //  derivative
          if (azz < 0) azz += 360.0;
          if (el < 0) el = 0;
          iaz = azz + 0.5;
          iel = el + 0.5;
          if (iaz < 0) iaz = 0;
          // if (iaz > 359) iaz = 359; SGM -- this is wrong, it should wrap around.
          if (iaz > 359) iaz -= 360;
          if (iel > 90) iel = 90;
          // frq150 = 1;  // choice makes no difference
          amp150 = sang * azel[(iaz + iel * 360) * NBEAM + frq150];
          mp = mapgal[i][j];

          //   mp = mp*1.1;  scale no effect
          //    if(site < 10) mp += -5;  // offset has a large effect -5 best
          //    low2
          //  Haslam 384 -site 0 Guz 384 site 10 Haslam+spind 256 site 100
          //  Guz+spind 256 site 1000
          aind = 2.5;
          if (el > 0 && !(skymode & 128)) {
            //    aind = log((map45[i][j]-cmbb)/(mp-cmbb))/log(408.0/45.0);
            aind = spind[i][j];
            //    if(aind > 3 || aind < 2) aind = 2.6;
            if (aind > 3) aind = 3.0;
            if (aind < 2) aind = 2.0;
          }
          //   aind = 2.6; // constant no effect
          for (frq = 0; frq < nbeam; frq++) {
            amp = sang * azel[(iaz + iel * 360) * NBEAM + frq];
            wsum = wsum2 = 0;
            if (el > 0) {
              fr = frqst + frq * frqspac;
              if (skymode & 128) {
                wsum = (mp - cmbb) * poww[frq] + cmbb;  // speeds up execution
                //                   if(fabs(glat)<10) wsum *=
                //                   pow(fr*0.5/408.0,0.1);  // correct galactic
                //                   plane spectral index difference
              } else {
                if (map >= 2) {
                  if (map == 2) wsum = (mp - cmbb) * pow(fr / 408.0, -aind) + cmbb;
                  if (map == 3) wsum = (mp - cmbb) * pow(fr / 45.0, -aind) + cmbb;
                } else {
                  gp = log(fr / 150.0);
                  if (fabs(glat) < 9.0)
                    wb = (cos(glon * PI / 180.0) + 1) / 2.0;
                  else
                    wb = 0;  // values of 10 deg and 2.4 not yet well determined
                  wsum = (mp - 3.0) * ((1 - wb) * pow(fr / 408.0, -2.52 + 0.016 * gp) + wb * pow(fr / 408.0, -2.8 - 0.13 * gp)) +
                         3.0;  // best values from 2015:024
                  if (fabs(glat - 45) < 15 && fabs(glon - 27) < 15)
                    wsum = mp * pow(150.0 / 408.0, -2.5) * pow(fr / 150.0,
                                                               -3.0 - 0.5 * gp);  // test north polar spur region
                }
              }
              if (skymode & 4) wsum += telec * opac * pow(fr / 150.0, -2.0);           // ion emission
              if (skymode & 8) wsum += (telec - wsum) * opac * pow(fr / 150.0, -2.0);  // ion emission and absorption
              if (skymode & 16) wsum = 500.0 * pow(fr / 150.0, -2.5);
              if (cmb & 16)
                wsum2 = 3e-3 * wsum *
                        cos(2 * PI * fr * 1e6 * 50.0 * (1.0 - cos(el * PI / 180.0) * cos((azz - 90) * PI / 180.0)) /
                            3e8);  // part in 1e3 scatter object east of antenna
              if (cmb & 256)
                wsum2 = -0.5 * (1 - exp(-7 * exp(-(fr - 78) * (fr - 78) * (-log(-log((1 + exp(-7.0)) / 2.0) / 7.0)) / (19 * 19 * 0.25)))) /
                        (1 - exp(-7.0));  // adding signature test
              //  printf("fr %f amp %e wsum %e wsum2 %e azz %lf el
              //  %lf\n",fr,amp,wsum,wsum2,azz,el);
            } else
              wsum = 300.0;
            if (wsum > max) max = wsum;
            if (frq==0){
              fprintf(file3, "%.15e %.15e %.15e %.15e %.15e %.15e %d %d %.15e\n", azz, el, wsum, amp, amp150, wsum2, iaz, iel, sang);
            }
            fsum[frq] += amp * (wsum + wsum2);
            fsum1[frq] += amp;
            fsum2[frq] += amp150 * wsum;
            
          }
        }
      }
      fclose(file3);

      sprintf(name, "beamfac%d.txt",n);
      file3 = fopen(name, "w");
      fprintf(file3, "frq bwfg bm bwfg_ref\n");
      for(frq=0;frq<nbeam;frq++){
        fprintf(file3, "%f %f %f\n", fsum[frq], fsum1[frq], fsum2[frq]);
      }
      fclose(file3);
    }
    //   if(nbeam > NFIT) bfit = NFIT;   // 30 needed when azelq spacing is 2
    //   MHz
    if (bfit >= 0) bfit = nbeam - bfit;  // input bfit = 1 useful to avoid ripple at ends of fit
    if (bfit < 0) bfit = -bfit;
    if (bfit > 50) bfit = 50;
    *bspac = (nbeam * 1.2) * frqspac;  // need +20% to avoid ripple
    for (i = 0; i < bfit; i++) {
      for (j = 0; j < nbeam; j++) {
        //      if(j*frqspac+frqst >= 100 && j*frqspac+frqst < 200)  wt[j] = 1;
        //      else wt[j] = 0;
        wt[j] = 1;
        if (bfit >= 7 && bfit <= 14)
          fitf[j + i * nbeam] = pow(((j * frqspac + frqst)) / (*freqref),
                                    i);  // increased to 12 but 9 normally used 12nov18
        else {
          if ((i % 2) == 0) fitf[j + i * nbeam] = cos((j * frqspac + frqst - (*freqref)) * 2 * PI * (i / 2) / (*bspac));
          if ((i % 2)) fitf[j + i * nbeam] = sin((j * frqspac + frqst - (*freqref)) * 2 * PI * ((i + 1) / 2) / (*bspac));
        }
      }
    }

    for (frq = 0; frq < nbeam; frq++) data[frq] = (fsum[frq] / fsum1[frq]) / (fsum2[frq] / fsum2[frq150]);  // normalize to 150 21 Jan 16
    if (cmb & 512)
      for (frq = 0; frq < nbeam; frq++)
        data[frq] = (fsum[frq] / fsum1[frq]) / pow((frqst + frq * frqspac) / ((double)(frqst + frq150 * frqspac)),
                                                   -2.5);  // previous
    polyfitr(bfit, nbeam, data, mcalc, wt, data2, fitf, aarr, bbrr);
    rms = 0;
    sum = 0;
    for (i = 0; i < nbeam; i++) {
      rms += wt[i] * (data[i] - data2[i]) * (data[i] - data2[i]);
      sum += wt[i];
    }
    //    for(i=0;i<nbeam;i++) printf("Spbeam data %d %f %f diff %e bbrr %e rms
    //    %f lst
    //    %5.2f\n",frqst+i*frqspac,data[i],data2[i],data[i]-data2[i],bbrr[i],sqrt(rms/sum),*lst);
    //    for(i=0;i<nbeam;i++) poww[i] = frqst+i*frqspac;
    //    plotfspec(nbeam,poww,data,data2,wt,0,0,6,name,name,name);
    printf("beamfit nbeam %d bfit %d rms %e max %e\n", nbeam, bfit, sqrt(rms / sum), max);
    for (i = 0; i < bfit; i++) bb[i] = bbrr[i];

    return bfit;
  }
  sum = 0;
  for (i = 0; i < bfit; i++) {
    if (bfit >= 7 && bfit <= 14)
      sum += bb[i] * pow(freq / (*freqref),
                         i);  // best since spectral index taken out cannot fit
                              // constant perfectly with pow(freq/150.0,-2.5+i);
    else {
      if ((i % 2) == 0) sum += bb[i] * cos((freq - (*freqref)) * 2 * PI * (i / 2) / (*bspac));
      if ((i % 2)) sum += bb[i] * sin((freq - (*freqref)) * 2 * PI * ((i + 1) / 2) / (*bspac));
    }
  }
  sum = sum * pow(freq / (*freqref), -2.5);
  return sum;
}

//  alternate for test

double beamcorr2(double freq, int bfit, double ang, double secs, double mcalc[], double fitf[], long double aarr[], double bbrr[], int skymode,
                 int integ, double *lst, double bb[], double *bspac, int cmb, double *freqref, double fbstart, double fbstop, int site) {
  char name[255], buf[32768], *p;
  FILE *file3;
  static double azel[360 * 91 * NBEAM], data[NBEAM], data2[NBEAM], wt[NBEAM], fsum[NBEAM], fsum1[NBEAM], fsum2[NBEAM];
  static float gmap[3145728];
  double k, sum, amp, amp150, raa, dec, lat, lon, gstt, glat, glon, wsum, azz, el, sang, mp, dp, gp, opac, telec, rms, wb;
  double ssecs, cmbb, theta, R, h, phi;
  int i, j, n, nn, az, frq, frq150, iaz, iel, m, frqst, frqspac, nbeam;
  double cost, sint, zb, antght, d;
  frqst = 0;
  frqspac = 0;
  n = 0;
  if (integ > 0) {
    if (cmb & 2)
      cmbb = 2.725;
    else
      cmbb = 0;
    i = -1;
    sprintf(name, "azelq.txt");
    if ((file3 = fopen(name, "r")) == NULL) {
      printf("%s error\n", name);
      return 0;
    }
    while (fgets(buf, 32768, file3) != 0) {
      if (strstr(buf, "#FREQ")) {
        sscanf(buf, "%*s %d", &frq);
        i++;
        if (i == 0 && fbstart < frq) frqst = frq;
        if (frq == fbstart) {
          frqst = frq;
          i = 0;
        }
        if (frq == fbstop) n = i;
        if (i == 1) frqspac = frq - frqst;
      }
      if (buf[0] == 'a' && i >= 0 && i < NBEAM) {
        p = &buf[3];
        j = 0;
        az = 0;
        az = strtod(p, &p) + ang;
        if (az < 0) az += 360;
        if (az < 0) az += 360;
        if (az >= 360) az -= 360;
        while (j < 91 && p && *p != '\n') {
          azel[(az + j * 360) * NBEAM + i] = pow(10.0, 0.1 * strtod(p, &p));
          //                           printf("az %d el %d fr %d amp
          //                           %f\n",az,j,frqst+i*frqspac,azel[(az+j*360)*NBEAM+i]);

          if (frqst < 70.0)
            *freqref = 75.0;
          else
            *freqref = 150.0;
          if (skymode >= 512) {
            el = j;
            cost = cos(el * PI / 180.0) * cos((az + ang) * PI / 180.0);
            zb = sin(el * PI / 180.0);
            sint = sqrt(1.0 - cost * cost);
            //  if(fabs(sint)>1e-6 && el>=0.0) amp=(cos(0.5*PI*cost))/sint;
            if (fabs(sint) > 1e-6 && el >= 0.0) {
              amp = PI * ((frqst + i * frqspac) / 150.0) / 2.0;  // half-wavelength at 150 MHz dipole
              amp = (cos(amp * cost) - cos(amp)) / sint;
            } else
              amp = 1e-6;
            antght = 0.25;  // height in wavelengths at 150 MHz
            if (zb > 0.0)
              amp = amp * sin(PI * 2.0 * zb * antght * (frqst + i * frqspac) / 150.0);
            else
              amp = 1e-6;
            // if((az==0 || az==90) && el==60)     printf("az %d el %d fr %d
            // amp1 %f amp2 %f ratio
            // %f\n",az,j,frqst+i*frqspac,azel[(az+j*360)*NBEAM+i],amp*amp,azel[(az+j*360)*NBEAM+i]/(amp*amp));
            azel[(az + j * 360) * NBEAM + i] = amp * amp;
          }
          j++;
        }
      }
    }
    fclose(file3);
    if (i < NBEAM + 1)
      nbeam = i + 1;
    else
      return 0;
    if (n) nbeam = n + 1;

    for (i = 0; i < 3145728; i++) gmap[i] = -1;
    if ((file3 = fopen("gmap.bin", "rb")) == NULL) {
      printf("gmap.bin not found\n");
      return 0;
    }
    if (!fread(gmap, sizeof(gmap), 1, file3)) printf("read error here %f %f\n", gmap[0], gmap[1]);
    fclose(file3);

    telec = 1e3;
    if (site == 0) {
      lon = 116.5 * PI / 180.0;  // Boolardy
      lat = -26.7 * PI / 180.0;  // Boolardy
    }
    for (m = 0; m < nbeam; m++) fsum[m] = fsum1[m] = fsum2[m] = 0;
    frq150 = 0;
    k = 1e6;
    for (frq = 0; frq < nbeam; frq++)
      if (fabs(frqst + frq * frqspac - (*freqref)) < k) {
        k = fabs(frqst + frq * frqspac - (*freqref));
        frq150 = frq;
      }  // find closest
    //    for(frq=0;frq<nbeam;frq++) if(fabs(frqst+frq*frqspac-150.0) < 0.1)
    //    frq150=frq;
    nn = integ / 1800 + 1;  // averaging over 3600 sec makes little difference to result
    d = (double)integ / (double)nn;
    for (n = 0; n < nn; n++) {
      ssecs = secs + d / 2.0 + n * d - ((double)integ) / 2.0;
      gstt = gst(ssecs);
      *lst = (gstt + lon) * 12.0 / PI;
      if (*lst < 0) *lst += 24.0;
      if (*lst > 24.0) *lst -= 24.0;
      for (i = 0; i < 3145728; i++) {
        pix2ang_ring(512, i, &theta, &phi);
        glon = phi * 180.0 / PI;
        glat = 90.0 - theta * 180.0 / PI;
        GalactictoRadec(glat, glon, &raa, &dec);
        radec_azel(gstt - raa + lon, dec, lat, &azz, &el);
        azz = azz * 180.0 / PI;
        el = el * 180.0 / PI;
        sang = (180.0 / 512.0) * cos(glat * PI / 180.0) * (360.0 / 1024.0) * PI * PI / (180.0 * 180.0);
        sang = 1;  // this should the case for HEALpix
        phi = el * PI / 180.0;
        R = 6356.0;
        h = 300.0;
        opac = (1.8 / 500.0) * (R + h) / sqrt(R * R * sin(phi) * sin(phi) + 2 * R * h + h * h);  // from derivative
        m = 0;
        if (azz < 0) azz += 360.0;
        if (el < 0) el = 0;
        iaz = azz + 0.5;
        iel = el + 0.5;
        if (iaz < 0) iaz = 0;
        if (iaz > 359) iaz = 359;
        if (iel > 90) iel = 90;
        amp150 = sang * azel[(iaz + iel * 360) * NBEAM + frq150];
        for (frq = 0; frq < nbeam; frq++) {
          amp = sang * azel[(iaz + iel * 360) * NBEAM + frq];
          if (skymode & 32) {
            if (iel > 5)
              amp = 1;
            else
              amp = 0;
          }  // isotropic antenna
          if (skymode & 64) {
            if (iel > 80)
              amp = 1;
            else
              amp = 0;
          }  // beam at Zenith
          if (el > 0) {
            mp = gmap[i];
            if (skymode & 1)
              dp = -0.05;
            else
              dp = 0;  // spectral index variation
            if (fabs(glat) > 20.0) dp = -dp;
            if (skymode & 2)
              gp = -0.12 * log((frqst + frq * frqspac) / 150.0);  // derivative Angelica's gamma
            else
              gp = 0.0;
            if (skymode & 128)
              wsum = (mp - cmbb) * pow((frqst + frq * frqspac) / 408.0, -2.5 + dp + gp) + cmbb;
            else {
              // correct 408 MHz map to 150 MHz - makes slight improvement at
              // most GHAs - see memo 7
              if (fabs(glat) < 10.0)
                wb = (cos(glon * PI / 180.0) + 1.0) / 2.0;
              else
                wb = 0;
              wsum = (mp - 3.0) * (wb * pow(150.0 / 408.0, -2.75 + 2.5) + (1 - wb)) * pow((frqst + frq * frqspac) / 408.0, -2.5 + dp + gp) + 3.0;
            }
            if (skymode & 4) wsum += telec * opac * pow((frqst + frq * frqspac) / 150.0, -2.0);  // ion emission
            if (skymode & 8) wsum += (telec - wsum) * opac * pow((frqst + frq * frqspac) / 150.0,
                                                                 -2.0);  // ion emission and absorption
            if (skymode & 16) wsum = 500.0 * pow((frqst + frq * frqspac) / 150.0, -2.5);
          } else
            wsum = 300.0;
          fsum[m] += amp * wsum;
          fsum1[m] += amp;
          fsum2[frq] += amp150 * wsum;
          m++;
        }
      }
    }

    //   if(nbeam > NFIT) bfit = NFIT;   // 30 needed when azelq spacing is 2
    //   MHz
    if (bfit >= 0) bfit = nbeam - bfit;  // input bfit = 1 useful to avoid ripple at ends of fit
    if (bfit < 0) bfit = -bfit;
    if (bfit > 50) bfit = 50;
    *bspac = (nbeam * 1.2) * frqspac;  // need +20% to avoid ripple
    for (i = 0; i < bfit; i++) {
      for (j = 0; j < nbeam; j++) {
        //      if(j*frqspac+frqst >= 100 && j*frqspac+frqst < 200)  wt[j] = 1;
        //      else wt[j] = 0;
        wt[j] = 1;
        if (bfit >= 7 && bfit <= 12)
          fitf[j + i * nbeam] = pow(((j * frqspac + frqst)) / (*freqref),
                                    i);  // increased to 12 but 9 normally used 12nov18
        else {
          if ((i % 2) == 0) fitf[j + i * nbeam] = cos((j * frqspac + frqst - (*freqref)) * 2 * PI * (i / 2) / (*bspac));
          if ((i % 2)) fitf[j + i * nbeam] = sin((j * frqspac + frqst - (*freqref)) * 2 * PI * ((i + 1) / 2) / (*bspac));
        }
      }
    }

    for (frq = 0; frq < nbeam; frq++) data[frq] = (fsum[frq] / fsum1[frq]) / (fsum2[frq] / fsum2[frq150]);  // normalize to 150 21 Jan 16
    //     for(frq=0;frq<nbeam;frq++)
    //     data[frq]=(fsum[frq]/fsum1[frq])/pow((frqst+frq*frqspac)/150.0,-2.5);
    //     // previous
    polyfitr(bfit, nbeam, data, mcalc, wt, data2, fitf, aarr, bbrr);
    rms = 0;
    sum = 0;
    for (i = 0; i < nbeam; i++) {
      rms += wt[i] * (data[i] - data2[i]) * (data[i] - data2[i]);
      sum += wt[i];
    }
    //    for(i=0;i<nbeam;i++) printf("Spbeam data %d %f %f diff %e bbrr %e rms
    //    %f lst
    //    %5.2f\n",frqst+i*frqspac,data[i],data2[i],data[i]-data2[i],bbrr[i],sqrt(rms/sum),*lst);
    //    for(i=0;i<nbeam;i++) poww[i] = frqst+i*frqspac;
    //    plotfspec(nbeam,poww,data,data2,wt,0,0,6,name,name,name);
    printf("beamfit nbeam %d bfit %d rms %e bbrr[0] %e fsum[0] %e\n", nbeam, bfit, sqrt(rms / sum), bbrr[0], fsum[0]);
    for (i = 0; i < bfit; i++) bb[i] = bbrr[i];

    return bfit;
  }
  sum = 0;
  for (i = 0; i < bfit; i++) {
    if (bfit >= 7 && bfit <= 12)
      sum += bb[i] * pow(freq / (*freqref),
                         i);  // best since spectral index taken out cannot fit
                              // constant perfectly with pow(freq/150.0,-2.5+i);
    else {
      if ((i % 2) == 0) sum += bb[i] * cos((freq - (*freqref)) * 2 * PI * (i / 2) / (*bspac));
      if ((i % 2)) sum += bb[i] * sin((freq - (*freqref)) * 2 * PI * ((i + 1) / 2) / (*bspac));
    }
  }
  sum = sum * pow(freq / (*freqref), -2.5);
  // printf("sum %f freq %f\n",sum,freq);
  return sum;
}
void GalactictoRadec(double glat, double glon, double *ra, double *dec)
/* galactic to radec  2000 epoch pole at 12h51.4 27.1 */
{
  double a, xg, yg, zg, xr, yr, zr, d0, dp, r0, rp;
  d0 = -(28.0 + 56.0 / 60.0) * PI / 180.0;
  r0 = (17.0 + 45.5 / 60.0) * PI / 12.0;
  dp = 27.1 * PI / 180.0;
  rp = (12.0 + 51.4 / 60.0) * PI / 12.0;
  zr = sin(d0);
  xr = cos(r0 - rp) * cos(d0);
  yr = sin(r0 - rp) * cos(d0);
  xg = xr * sin(dp) - zr * cos(dp);
  yg = yr;
  a = atan2(yg, xg);
  xg = cos((glon * PI / 180.0) + a) * cos(glat * PI / 180.0);
  yg = sin((glon * PI / 180.0) + a) * cos(glat * PI / 180.0);
  zg = sin(glat * PI / 180.0);
  xr = xg * sin(dp) + zg * cos(dp);
  yr = yg;
  zr = zg * sin(dp) - xg * cos(dp);
  *dec = atan2(zr, sqrt(xr * xr + yr * yr));
  *ra = atan2(yr, xr) + rp;
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

double polyfitrc(int npoly, int nfreq, double ddata[], double mcalc[], double wtt[], double dataout[], double fitfn[], long double aarr[],
                 double bbrr[], int mode) {
  int i, iter;
  double cons[50];
  polyfitr(npoly, nfreq, ddata, mcalc, wtt, dataout, fitfn, aarr, bbrr);
  if (mode == 0)
    return bbrr[0];
  else {
    for (iter = 0; iter < 10; iter++) {
      for (i = 0; i < npoly; i++)
        if (bbrr[i] < 0)
          cons[i] = mode * 1e-3;
        else
          cons[i] = 0;
      polyfitr2(npoly, nfreq, ddata, mcalc, wtt, dataout, fitfn, aarr, bbrr, cons);
    }
  }
  return bbrr[0];
}

double polyfitr2(int npoly, int nfreq, double ddata[], double mcalc[], double wtt[], double dataout[], double fitfn[], long double aarr[],
                 double bbrr[], double cons[]) {
  int i, j, k, kk, m1, m2;
  double re;
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
      if (wtt[k] > 0) re += mcalc[m1] * ddata[k] * wtt[k];
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
  for (i = 0; i < npoly; i++) aarr[i + i * npoly] += cons[i];  // add constraints
  qrd(aarr, npoly, bbrr);
  for (i = 0; i < nfreq; i++) {
    re = 0.0;
    for (j = 0; j < npoly; j++) { re += bbrr[j] * fitfun(j, i, nfreq, fitfn); }
    dataout[i] = re;
  }
  return bbrr[0];
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

      //     if(i==j) re += 1e-20;     // matrix inversion problem
      aarr[k] = re;
      //  printf("k %d %d %d aarr %Le\n",k,j,i,aarr[k]);
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

double fitfun(int j, int i, int nfreq, double fitfn[]) { return fitfn[i + j * nfreq]; }

void plotfspec(int np, double freqq[], double data[], double data2[], double wtt[], double ddscale, int mode, int num, char title[], char datfile[],
               char info[])
// plot the spectrum
{
  char txt[256];
  int k, iter, nter;
  double x, y, xp, yp, dmax, dmin, f, totpp, scale, step, h, s, b;
  double xoffset, yoffset, fstart, fstop, dscale, rms, rms1, av, avv, ss;
  FILE *file;

  sprintf(txt, "spe%d.pos", num);
  if ((file = fopen(txt, "w")) == NULL) {
    printf("cannot open %s:\n", txt);
    return;
  }
  fprintf(file, "%%!PS-Adobe-\n%c%cBoundingBox:  0 0 612 700\n%c%cEndProlog\n", '%', '%', '%', '%');
  fprintf(file, "1 setlinewidth\n");
  fprintf(file, "/Times-Roman findfont\n 14 scalefont\n setfont\n");

  xoffset = 80.0;
  yoffset = 100.0;
  fstart = freqq[0];
  fstop = freqq[np - 1];
  dmax = -1e99;
  dmin = 1e99;
  dscale = -1e99;
  for (k = 0; k < np; k++)
    if (wtt[k] && data2[k] > dmax) dmax = data2[k];
  for (k = 0; k < np; k++)
    if (wtt[k] && data2[k] < dmin) dmin = data2[k];
  for (k = 0; k < np; k++)
    if (wtt[k] && fabs(data2[k] - data[k]) > dscale) dscale = fabs(data2[k] - data[k]);
  if (dscale < 1e-6) dscale = 1;
  if (ddscale) dscale = ddscale;
  rms = 0;
  ss = 1e-30;
  av = avv = 0;
  for (k = 0; k < np; k++)
    if (wtt[k]) {
      rms += (data2[k] - data[k]) * (data2[k] - data[k]);
      av += data2[k] - data[k];
      avv += data2[k];
      ss++;
    }
  rms = sqrt(rms / ss);
  av = av / ss;
  avv = avv / ss;
  rms1 = 0;
  for (k = 0; k < np; k++)
    if (wtt[k]) { rms1 += (data2[k] - data[k] - av) * (data2[k] - data[k] - av); }
  rms1 = sqrt(rms1 / ss);
  nter = 3;
  dmin = 0;
  dmax = dmax * 1.4;
  k = dmax / 1000.0 + 0.5;
  dmax = k * 1000.0;
  //    if(dmax < 1000.0) dmax = 1000.0;
  //    if(dmax > 2000) dmax=2000;   // temporary fix
  if (dmax < 400.0) dmax = 400.0;
  if (mode) {
    dmax = 120;
    dmin = 0;
  }
  if (num == 5) {
    dmax = 5;
    for (k = 0; k < np; k++)
      if (data2[k] > dmax) dmax = data2[k];
    if (dmax > 10) dmax = 10;
  }
  scale = dmax - dmin;
  //    printf("dmax %f dmin %f\n",dmax,dmin);
  for (y = 0; y < 2; y++)
    fprintf(file,
            "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
            "sethsbcolor stroke\n",
            xoffset, y * 480 + yoffset, xoffset + 400.0, y * 480 + yoffset);
  fprintf(file,
          "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
          "sethsbcolor stroke\n",
          xoffset, 400 + yoffset, xoffset + 400.0, 400 + yoffset);
  fprintf(file,
          "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
          "sethsbcolor stroke\n",
          xoffset, yoffset, xoffset, 480.0 + yoffset);
  fprintf(file,
          "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
          "sethsbcolor stroke\n",
          xoffset + 400.0, yoffset, xoffset + 400.0, 480.0 + yoffset);
  for (iter = 0; iter < nter; iter++) {
    if (iter < 2) fprintf(file, "%d setlinewidth\n", 0);
    if (iter == 2) fprintf(file, "%d setlinewidth\n", 1);
    yp = 0;
    xp = 0;
    totpp = 0;
    for (k = 0; k < np; k++) {
      h = s = b = 0;
      x = (freqq[k] - fstart) * 400.0 / (fstop - fstart);
      totpp = (data[k] - dmin) / scale;
      if (iter == 1) totpp = (data2[k] - dmin) / scale;
      y = totpp * 400.0;
      if (y > 400.0) y = 400.0;
      if (y < 0.0) y = 0.0;
      if (iter == 2) {
        totpp = (data[k] - data2[k]) / dscale;
        y = totpp * 40.0 + 440.0;
        if (y > 480.0) y = 480;
        if (y < 400.0) y = 400;
      }
      if (k == 0 || x < xp) {
        xp = x;
        yp = y;
      }
      if (!iter || (iter == 1 && data2[k] < dmax) || iter == 2) {
        if (iter == 1 && k > 0 && wtt[k] == 0)
          fprintf(file, "%d setlinewidth\n", 1);
        else if (iter == 1)
          fprintf(file, "%d setlinewidth\n", 2);
        h = s = b = 0;
        if (iter < 2 || (iter == 2 && wtt[k] && wtt[k - 1] && wtt[k + 1]))
          fprintf(file,
                  "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n %5.2f "
                  "%5.2f %5.2f sethsbcolor stroke\n",
                  xp + xoffset, yp + yoffset, x + xoffset, y + yoffset, h, s, b);
        xp = x;
        yp = y;
      }
    }
  }
  fprintf(file, "1 setlinewidth\n");
  step = 5.0;
  if (fstop - fstart > 60) step = 10;
  if (fstop - fstart > 100) step = 20;
  for (f = (int)fstart; f <= fstop + step * 0.1; f += step) {
    x = (f - fstart) * 400.0 / (fstop - fstart);
    y = 0;
    fprintf(file,
            "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
            "sethsbcolor stroke\n",
            x + xoffset, y + yoffset, x + xoffset, y + 5.0 + yoffset);
    sprintf(txt, "%5.1f", f);
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 15.0, y - 15.0 + yoffset, txt);
  }
  {
    for (f = dmin; f < dmax * 0.95; f += (dmax - dmin) / 10.0) {
      x = 0;
      y = (f - dmin) * 400.0 / (dmax - dmin);
      fprintf(file,
              "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
              "sethsbcolor stroke\n",
              x + xoffset, y + yoffset, x + xoffset + 5, y + yoffset);
      sprintf(txt, "%4.0f K", f);
      if (dmax < 10.0) sprintf(txt, "%3.1f K", f);
      fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 45.0, y - 1.0 + yoffset, txt);
    }
    for (f = -dscale; f <= dscale; f += dscale * 0.5) {
      x = 0;
      y = f * 40.0 / dscale + 440.0;
      fprintf(file,
              "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
              "sethsbcolor stroke\n",
              x + xoffset, y + yoffset, x + xoffset + 5, y + yoffset);
      if (dscale >= 10.0) sprintf(txt, "%4.0f K", f);
      if (dscale < 10.0) sprintf(txt, "%3.1f K", f);
      if (dscale < 1.0) sprintf(txt, "%4.0fmK", f * 1e3);
      fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 50.0, y - 1.0 + yoffset, txt);
    }
  }
  sprintf(txt, " Temperature (K)");
  fprintf(file, "%6.2f %6.2f moveto\n 90 rotate\n (%s) show\n -90 rotate\n", xoffset - 52.0, 265.0, txt);

  sprintf(txt, "Frequency (MHz)");
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 160.0, 70.0, txt);
  if (num == 0)
    sprintf(txt,
            "raw sky spectrum (thin line) vs calibrated spectrum of sky (thick "
            "line)");
  if (num == 1)
    sprintf(txt,
            "hot temperature (thin line) vs calibrated spectrum of hot "
            "load(thick line)");
  if (num == 2)
    sprintf(txt,
            "ambient temperature (thin line) vs calibrated spectrum of ambient "
            "load(thick line)");
  if (num == 3) sprintf(txt, "cable spectrum (thin line) vs fit to cable spectrum (thick line)");
  if (num == 4)
    sprintf(txt,
            "calibrated sky spectrum (thin line) vs sky model spectrum (thick "
            "line)");
  if (num == 5) sprintf(txt, "spectrum of loss (thin line) vs spectrum of loss (thick line)");
  if (num == 6) sprintf(txt, "spectrum of beam (thin line) vs spectrum of fit (thick line)");
  if (num == 99)
    sprintf(txt,
            "calibrated sky spectrum (thin line) vs fit to sky spectrum (thick "
            "line) %s",
            info);
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset, 55.0, txt);
  sprintf(txt, "%s rms diff. %5.3f %5.3f av %5.0f K %s", title, rms, rms1, avv, datfile);
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset - 10.0, 40.0, txt);
  fprintf(file, "showpage\n%c%cTrailer\n", '%', '%');
  fclose(file);
}

void outfspec(int np, double freqq[], double data[], double skymodel[], double sunn[], double lst, int mode, int nline, char title[])
// writeout calibrated spectrum
{
  char txt[256];
  int k;
  FILE *file;

  sprintf(txt, "spe0.txt");
  if ((file = fopen(txt, "w")) == NULL) {
    printf("cannot open %s:\n", txt);
    return;
  }

  for (k = 0; k < np; k++) {
    if (mode == 1) {
      fprintf(file,
              "freq %10.6f tantenna %10.6f K skymodel %10.6f K ant_area %8.6f "
              "m^2 Tsun_per_SFU %5.3f K lst %5.2f %s\n",
              freqq[k], data[k], skymodel[k], sunn[k], sunn[k] * 1e-22 / (2.0 * 1.38e-23), lst, title);
    }
    if (mode == 0) {
      fprintf(file, "freq %12.6f tantenna %12.6f K skymodel %12.6f K wt %1.0f %s %5d\n", freqq[k], data[k], skymodel[k], sunn[k], title, nline);
    }
  }
  fclose(file);
}

void outcal(int np, double freqq[], complex double s11lna[], double sca[], double ofs[], double tlnau[], double tlnacc[], double tlnacs[],
            double wtcal[], char title[])
// writeout calibrated spectrum
{
  char txt[256];
  int k;
  FILE *file;

  sprintf(txt, "specal.txt");
  if ((file = fopen(txt, "w")) == NULL) {
    printf("cannot open %s:\n", txt);
    return;
  }

  for (k = 0; k < np; k++) {
    fprintf(file,
            "freq %1.16e s11lna %1.16e %1.16e sca %1.16e ofs %1.16e tlnau "
            "%1.16e tlnac %1.16e tlnas %1.16e wtcal %1.0f %s\n",
            freqq[k], creal(s11lna[k]), cimag(s11lna[k]), sca[k], ofs[k], tlnau[k], tlnacc[k], tlnacs[k], wtcal[k], title);
  }
  fclose(file);
}

void outsim(int np, double freqq[], double spant[], int yr, int dy, int hr, int mn, int sc)
// simulated spectrum
{
  char txt[256];
  int k;
  FILE *file;

  sprintf(txt, "spesim.txt");
  if ((file = fopen(txt, "w")) == NULL) {
    printf("cannot open %s:\n", txt);
    return;
  }

  for (k = 0; k < np; k++) {
    if (k == 0)
      fprintf(file, "%10.6f %10.6f 1 3 // sim %04d:%03d:%02d:%02d:%02d 0  0 0 0 0\n", freqq[k], spant[k], yr, dy, hr, mn, sc);
    else
      fprintf(file, "%10.6f %10.6f 1\n", freqq[k], spant[k]);
  }
  fclose(file);
}

void readcal(int *np, double freqq[], complex double s11lna[], double sca[], double ofs[], double tlnau[], double tlnacc[], double tlnacs[],
             double wtcal[]) {
  char txt[256], buf[32768];
  double freq, re, im, ssca, oofs, t0, t1, t2, t3;
  complex double T;
  int k;
  FILE *file;

  *np = 0;
  sprintf(txt, "specal.txt");
  if ((file = fopen(txt, "r")) == NULL) {
    printf("cannot open %s:\n", txt);
    return;
  }
  k = 0;
  while (fgets(buf, 32768, file) != 0) {
    sscanf(buf, "%*s %lf %*s %lf %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf", &freq, &re, &im, &ssca, &oofs, &t0, &t1, &t2, &t3);
    T = re + im * I;
    //                    T =
    //                    T*cexp(-2.0*PI*freq*1e6*(0.0e-12)*I)*pow(10.0,0.05*0.2);
    //                    // for test corrections
    freqq[k] = freq;
    s11lna[k] = T;
    sca[k] = ssca;
    ofs[k] = oofs;
    tlnau[k] = t0;
    tlnacc[k] = t1;
    tlnacs[k] = t2;
    wtcal[k] = t3;
    k++;
  }
  *np = k;
  for (k = 0; k < *np; k++) {
    //    printf("freq %10.6f s11lna %10.6f %10.6f sca %10.6f ofs %10.6f tlnau
    //    %10.6f tlnac %10.6f tlnas %10.6f wtcal %1.0f\n",
    //       freqq[k],creal(s11lna[k]),cimag(s11lna[k]),sca[k],ofs[k],tlnau[k],tlnacc[k],tlnacs[k],wtcal[k]);
  }
  fclose(file);
}

void plotvna(int np, double freqq[], complex double data[], double wtt[], int npout, double freqq2[], complex double data2[], double dbdiff[],
             double pdiff[], int mode, double rmsdb, double rmsphase, int nfit, char *fname)
// plot the spectrum
{
  char txt[256];
  int k, iter, nn, nfirst;
  double x, y, xp, yp, dmax, dmin, f, totpp, scale, step, h, s, b, dscale, pscale;
  double xoffset, yoffset, fstart, fstop;  // ypos[4];
  FILE *file;

  sprintf(txt, "spevna%d.pos", mode);
  if ((file = fopen(txt, "w")) == NULL) {
    printf("cannot open %s:\n", txt);
    return;
  }
  fprintf(file, "%%!PS-Adobe-\n%c%cBoundingBox:  0 0 612 700\n%c%cEndProlog\n", '%', '%', '%', '%');
  fprintf(file, "1 setlinewidth\n");
  fprintf(file, "/Times-Roman findfont\n 14 scalefont\n setfont\n");

  xoffset = 80.0;
  yoffset = 100.0;
  fstart = freqq2[0];
  fstop = freqq2[npout - 1];
  dmax = -1e99;
  dmin = 1e99;
  for (k = 0; k < npout; k++)
    if (10.0 * log10(creal(data2[k] * conj(data2[k]))) > dmax) dmax = 10.0 * log10(creal(data2[k] * conj(data2[k])));
  for (k = 0; k < npout; k++)
    if (10.0 * log10(creal(data2[k] * conj(data2[k]))) < dmin) dmin = 10.0 * log10(creal(data2[k] * conj(data2[k])));
  dscale = pscale = 1e-99;
  for (k = 0; k < np; k++)
    if (fabs(dbdiff[k]) > dscale) dscale = fabs(dbdiff[k]);
  for (k = 0; k < np; k++)
    if (fabs(pdiff[k]) > pscale) pscale = fabs(pdiff[k]);
  k = dmin * 0.1 - 1.0;
  dmin = k * 10;
  k = dmax * 0.1 + 0.5;
  dmax = k * 10;
  scale = dmax - dmin;
  if (scale == 0.0) {
    dmax = 1;
    dmin = -1;
    scale = 2.0;
  }
  //   printf("dmax %f dmin %f np %d npout %d\n",dmax,dmin,np,npout);
  for (y = 0; y < 2; y++)
    fprintf(file,
            "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
            "sethsbcolor stroke\n",
            xoffset, y * 480 + yoffset, xoffset + 400.0, y * 480 + yoffset);
  fprintf(file,
          "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
          "sethsbcolor stroke\n",
          xoffset, 400 + yoffset, xoffset + 400.0, 400 + yoffset);
  fprintf(file,
          "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
          "sethsbcolor stroke\n",
          xoffset, yoffset, xoffset, 480.0 + yoffset);
  fprintf(file,
          "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
          "sethsbcolor stroke\n",
          xoffset + 400.0, yoffset, xoffset + 400.0, 480.0 + yoffset);
  fprintf(file,
          "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
          "sethsbcolor stroke\n",
          xoffset + 0.0, yoffset, xoffset + 0.0, 480.0 + yoffset);
  for (iter = 0; iter < 6; iter++) {
    if (iter % 2)
      fprintf(file, "%d setlinewidth\n", 1);  // for phase
    else
      fprintf(file, "%d setlinewidth\n", 0);
    if (iter < 2 || iter > 3) fprintf(file, "%d setlinewidth\n", 1);
    yp = 0;
    xp = 0;
    totpp = 0;
    if (iter < 2 || iter > 3)
      nn = np;
    else
      nn = npout;
    nfirst = 1;
    for (k = 0; k < nn; k++) {
      h = s = b = 0;
      if (iter < 2 || iter > 3)
        x = (freqq[k] - fstart) * 400.0 / (fstop - fstart);
      else
        x = (freqq2[k] - fstart) * 400.0 / (fstop - fstart);
      if (iter == 0) totpp = (10.0 * log10(creal(data[k] * conj(data[k]))) - dmin) / scale;
      if (iter == 2) totpp = (10.0 * log10(creal(data2[k] * conj(data2[k]))) - dmin) / scale;
      if (iter == 1) totpp = (atan2(cimag(data[k]), creal(data[k])) + PI) / (2.0 * PI);
      if (iter == 3) totpp = (atan2(cimag(data2[k]), creal(data2[k])) + PI) / (2.0 * PI);
      if (iter == 4) totpp = dbdiff[k] / dscale;
      if (iter == 5) totpp = pdiff[k] / pscale;
      if (iter < 4) {
        y = totpp * 400.0;
        if (y > 400.0) y = 400.0;
      }
      if (iter == 4) y = totpp * 8.0 + 450;
      if (iter == 5) y = totpp * 8.0 + 430;
      if (x >= -0.5 && x <= 400.5 && (iter < 4 || (wtt[k] && iter > 3))) {
        if (nfirst) {
          nfirst = 0;
          xp = x;
          yp = y;
        }
        //             if(fabs(k-0.7*nn) < 1.0) ypos[iter]=y+10.0;
        if (iter < 2)
          fprintf(file, "newpath\n %5.1f %5.1f %5.1f 0 360 arc\n closepath\n stroke\n", x + xoffset, y + yoffset, 1.0 + iter);
        else
          fprintf(file,
                  "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n %5.2f "
                  "%5.2f %5.2f sethsbcolor stroke\n",
                  xp + xoffset, yp + yoffset, x + xoffset, y + yoffset, h, s, b);
        xp = x;
        yp = y;
      }
    }
  }
  fprintf(file, "1 setlinewidth\n");
  step = 5.0;
  if (fstop - fstart > 60.0) step = 10.0;
  if (fstop - fstart > 100.0) step = 20.0;
  for (f = (int)fstart; f <= fstop + step * 0.1; f += step) {
    x = (f - fstart) * 400.0 / (fstop - fstart);
    y = 0.0;
    fprintf(file,
            "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
            "sethsbcolor stroke\n",
            x + xoffset, y + yoffset, x + xoffset, y + 5.0 + yoffset);
    sprintf(txt, "%5.1f", f);
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 15.0, y - 15.0 + yoffset, txt);
  }
  for (f = dmin; f <= dmax; f += (dmax - dmin) / 10.0) {
    x = 0;
    y = (f - dmin) * 400.0 / (dmax - dmin);
    fprintf(file,
            "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
            "sethsbcolor stroke\n",
            x + xoffset, y + yoffset, x + xoffset + 5, y + yoffset);
    {
      sprintf(txt, "%4.0f", f);
      fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 25.0, y - 1.0 + yoffset, txt);
    }
  }
  for (f = -180.0; f <= 180.0; f += 30.0) {
    x = 0;
    y = (f + 180.0) * 400.0 / 360.0;
    fprintf(file,
            "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
            "sethsbcolor stroke\n",
            x + xoffset + 400.0, y + yoffset, x + xoffset + 395.0, y + yoffset);
    {
      sprintf(txt, "%4.0f", f);
      fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset + 405.0, y - 1.0 + yoffset, txt);
    }
  }
  {
    sprintf(txt, "fit %5.4f p-p dB", 2.0 * dscale);
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 50.0, yoffset + 465.0, txt);
  }
  {
    sprintf(txt, "fit %5.4f p-p deg", 2.0 * pscale);
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 50.0, yoffset + 405.0, txt);
  }
  {
    sprintf(txt, "magnitude(thin line)");
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 280.0, yoffset + 25.0, txt);
  }
  {
    sprintf(txt, "phase(thick line)");
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 280.0, yoffset + 10.0, txt);
  }
  sprintf(txt, " S11 magnitude (dB)");
  fprintf(file, "%6.2f %6.2f moveto\n 90 rotate\n (%s) show\n -90 rotate\n", xoffset - 35.0, 265.0, txt);
  sprintf(txt, " S11 phase (degrees)");
  fprintf(file, "%6.2f %6.2f moveto\n 90 rotate\n (%s) show\n -90 rotate\n", xoffset + 440.0, 265.0, txt);
  sprintf(txt, "Frequency (MHz)");
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 180.0, 65.0, txt);
  if (mode == 0) sprintf(txt, "%d term Fit to antenna S11", nfit);
  if (mode == 1) sprintf(txt, "%d term Fit to LNA  S11", nfit);
  if (mode == 2) sprintf(txt, "%d term Fit to hot S11", nfit);
  if (mode == 3) sprintf(txt, "%d term Fit to amb S11", nfit);
  if (mode == 4) sprintf(txt, "%d term Fit to open cab S11", nfit);
  if (mode == 5) sprintf(txt, "%d term Fit to shorted cab S11", nfit);
  fprintf(file, "/Times-Roman findfont\n 18 scalefont\n setfont\n");
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 30.0, 45.0, txt);
  fprintf(file, "/Times-Roman findfont\n 14 scalefont\n setfont\n");
  sprintf(txt, "rms diff %6.3f dB %6.3f deg\n", rmsdb, rmsphase);
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 250.0, 45.0, txt);
  fprintf(file, "%6.2f %6.2f moveto\n (file: %s) show\n", xoffset + 50.0, 30.0, fname);
  fprintf(file, "showpage\n%c%cTrailer\n", '%', '%');
  fclose(file);
}

void plotnwave(int np, double freqq[], double data[], double data2[], double data3[], double delaycab, int wfit)
// plot the spectrum
{
  char txt[256];
  int k, iter, nter;
  double x, y, xp, yp, dmax, dmin, f, totpp, scale, step, h, s, b, dly, cc, cs;
  double xoffset, yoffset, fstart, fstop, ypos[3];
  FILE *file;

  if ((file = fopen("spewav.pos", "w")) == NULL) {
    printf("cannot open spewav.pos:\n");
    return;
  }
  fprintf(file, "%%!PS-Adobe-\n%c%cBoundingBox:  0 0 612 700\n%c%cEndProlog\n", '%', '%', '%', '%');
  fprintf(file, "1 setlinewidth\n");
  fprintf(file, "/Times-Roman findfont\n 14 scalefont\n setfont\n");

  xoffset = 80.0;
  yoffset = 100.0;
  fstart = freqq[0];
  fstop = freqq[np - 1];
  dmax = 100;
  dmin = 0;
  if (fstart < 70) dmax = 200;  // lowband
  nter = 3;
  scale = dmax - dmin;
  //    printf("dmax %f dmin %f\n",dmax,dmin);
  for (y = 0; y < 2; y++)
    fprintf(file,
            "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
            "sethsbcolor stroke\n",
            xoffset, y * 480 + yoffset, xoffset + 400.0, y * 480 + yoffset);
  fprintf(file,
          "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
          "sethsbcolor stroke\n",
          xoffset, yoffset, xoffset, 480.0 + yoffset);
  fprintf(file,
          "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
          "sethsbcolor stroke\n",
          xoffset + 400.0, yoffset, xoffset + 400.0, 480.0 + yoffset);
  fprintf(file,
          "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
          "sethsbcolor stroke\n",
          xoffset + 0.0, yoffset, xoffset + 0.0, 480.0 + yoffset);
  for (iter = 0; iter < nter; iter++) {
    if (iter < 2)
      fprintf(file, "%d setlinewidth\n", 1);
    else
      fprintf(file, "%d setlinewidth\n", 2);
    yp = 0;
    xp = 0;
    totpp = 0;
    for (k = 0; k < np; k++) {
      h = s = b = 0;
      x = (freqq[k] - fstart) * 400.0 / (fstop - fstart);
      totpp = (data[k] - dmin) / scale;
      if (iter == 1) totpp = (sqrt(data2[k] * data2[k] + data3[k] * data3[k]) - dmin) / scale;
      if (iter == 2) {
        dly = delaycab * 1e6;
        cc = data2[k] * cos(2.0 * PI * freqq[k] * dly) + data3[k] * sin(2.0 * PI * freqq[k] * dly);
        cs = data3[k] * cos(2.0 * PI * freqq[k] * dly) - data2[k] * sin(2.0 * PI * freqq[k] * dly);
        totpp = (atan2(cs, cc) + PI) / (2.0 * PI);
      }
      y = totpp * 480.0;
      if (y > 480.0) y = 480.0;
      if (k == 0) {
        xp = x;
        yp = y;
        ypos[iter] = y + 10;
      }
      fprintf(file,
              "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n %5.2f %5.2f "
              "%5.2f sethsbcolor stroke\n",
              xp + xoffset, yp + yoffset, x + xoffset, y + yoffset, h, s, b);
      xp = x;
      yp = y;
    }
    if (ypos[iter] > 470) ypos[iter] = 470;
  }
  fprintf(file, "1 setlinewidth\n");
  step = 5.0;
  if (fstop - fstart > 60.0) step = 10.0;
  if (fstop - fstart > 100.0) step = 20.0;
  for (f = (int)fstart; f <= fstop + step * 0.1; f += step) {
    x = (f - fstart) * 400.0 / (fstop - fstart);
    y = 0.0;
    fprintf(file,
            "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
            "sethsbcolor stroke\n",
            x + xoffset, y + yoffset, x + xoffset, y + 5.0 + yoffset);
    sprintf(txt, "%5.1f", f);
    fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 15.0, y - 20.0 + yoffset, txt);
  }
  for (f = dmin; f < dmax; f += (dmax - dmin) / 10.0) {
    x = 0;
    y = (f - dmin) * 480.0 / (dmax - dmin);
    fprintf(file,
            "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
            "sethsbcolor stroke\n",
            x + xoffset, y + yoffset, x + xoffset + 5, y + yoffset);
    {
      sprintf(txt, "%4.0f", f);
      fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset - 25.0, y - 1.0 + yoffset, txt);
    }
  }
  for (f = -180.0; f <= 180.0; f += 30.0) {
    x = 0;
    y = (f + 180.0) * 480.0 / 360.0;
    fprintf(file,
            "newpath\n %6.2f %6.2f moveto\n %6.2f %6.2f lineto\n 0 0 0 "
            "sethsbcolor stroke\n",
            x + xoffset + 400.0, y + yoffset, x + xoffset + 395.0, y + yoffset);
    {
      sprintf(txt, "%4.0f", f);
      fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", x + xoffset + 405.0, y - 1.0 + yoffset, txt);
    }
  }
  sprintf(txt, "magnitude of correlated wave");
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 190.0, yoffset + ypos[1], txt);
  sprintf(txt, "amplitude of uncorrelated wave");
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 190.0, yoffset + ypos[0], txt);
  sprintf(txt, "phase of correlated wave (thick line)");
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 190.0, yoffset + ypos[2], txt);
  sprintf(txt, " noise wave amplitude (K)");
  fprintf(file, "%6.2f %6.2f moveto\n 90 rotate\n (%s) show\n -90 rotate\n", xoffset - 35.0, 265.0, txt);
  sprintf(txt, " phase of correlated portion (degrees)");
  fprintf(file, "%6.2f %6.2f moveto\n 90 rotate\n (%s) show\n -90 rotate\n", xoffset + 440.0, 265.0, txt);
  sprintf(txt, "Frequency (MHz)");
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 160.0, 65.0, txt);
  sprintf(txt, "%5.5f ns delay taken out %d term Fit to noise waves", delaycab * 1e9, wfit);
  fprintf(file, "%6.2f %6.2f moveto\n (%s) show\n", xoffset + 20.0, 50.0, txt);

  fprintf(file, "showpage\n%c%cTrailer\n", '%', '%');
  fclose(file);
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

void radec_azel2(double ha, double sindec, double cosdec, double sinlat, double coslat, double *azs, double *elevs) {
  double w, r, zen, north;
  w = sin(ha) * cosdec;
  r = cos(ha) * cosdec;
  zen = r * coslat + sindec * sinlat;
  north = -r * sinlat + sindec * coslat;
  *elevs = atan2(zen, sqrt(north * north + w * w));
  *azs = atan2(-w, north);
  if (*azs < 0) *azs = *azs + TWOPI;
}

void pix2ang_ring(const long nside, long ipix, double *theta, double *phi) {
  /*
    from HEALPix  c-code
    c=======================================================================
    c     gives theta and phi corresponding to pixel ipix (RING)
    c     for a parameter nside
    c=======================================================================
  */

  int nl2, nl4, npix, ncap, iring, iphi, ip, ipix1;
  double fact1, fact2, fodd, hip, fihip;

  int ns_max = 8192;

  if (nside < 1 || nside > ns_max) {
    fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
    exit(0);
  }
  npix = 12 * nside * nside;  // ! total number of points
  if (ipix < 0 || ipix > npix - 1) {
    fprintf(stderr, "%s (%d): ipix out of range: %ld\n", __FILE__, __LINE__, ipix);
    exit(0);
  }

  ipix1 = ipix + 1;  // in {1, npix}
  nl2 = 2 * nside;
  nl4 = 4 * nside;
  ncap = 2 * nside * (nside - 1);  // ! points in each polar cap, =0 for nside =1
  fact1 = 1.5 * nside;
  fact2 = 3.0 * nside * nside;

  if (ipix1 <= ncap) {  //! North Polar cap -------------

    hip = ipix1 / 2.;
    fihip = floor(hip);
    iring = (int)floor(sqrt(hip - sqrt(fihip))) + 1;  // ! counted from North pole
    iphi = ipix1 - 2 * iring * (iring - 1);

    *theta = acos(1. - iring * iring / fact2);
    *phi = (1. * iphi - 0.5) * PI / (2. * iring);
  } else if (ipix1 <= nl2 * (5 * nside + 1)) {  // then ! Equatorial region ------

    ip = ipix1 - ncap - 1;
    iring = (int)floor(ip / nl4) + nside;  // ! counted from North pole
    iphi = (int)fmod(ip, nl4) + 1;

    fodd = 0.5 * (1 + fmod((double)(iring + nside),
                           2));  //  ! 1 if iring+nside is odd, 1/2 otherwise
    *theta = acos((nl2 - iring) / fact1);
    *phi = (1. * iphi - fodd) * PI / (2. * nside);
  } else {  //! South Polar cap -----------------------------------

    ip = npix - ipix1 + 1;
    hip = ip / 2.;
    /* bug corrige floor instead of 1.* */
    fihip = floor(hip);
    iring = (int)floor(sqrt(hip - sqrt(fihip))) + 1;  //     ! counted from South pole
    iphi = (int)(4. * iring + 1 - (ip - 2. * iring * (iring - 1)));

    *theta = acos(-1. + iring * iring / fact2);
    *phi = (1. * iphi - 0.5) * PI / (2. * iring);
  }
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
