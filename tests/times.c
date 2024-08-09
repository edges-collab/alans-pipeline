#include <math.h>
#include <stdio.h>

#define PI 3.1415926536
#define TWOPI 6.28318530717958

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

int main(){
    int i;
    int j;
    FILE *fp;
    fp = fopen("coords.txt", "w");
    fprintf(fp, "# LAT LON RA DEC AZ EL\n");

    double ssecs=1474043217.3333333;
    double sunra,sundec;
    double sunaz, sunel;
    double azz, el, raa, dec;

    double lon = 116.5 * PI / 180.0;  // Boolardy
    double lat = -26.7 * PI / 180.0;  // EDGES -26.72 116.61

    double gstt = gst(ssecs);
    sunradec(ssecs, &sunra, &sundec);
    radec_azel(gstt - sunra + lon, sundec, lat, &sunaz, &sunel);

    for(i=0;i<512;i++){
        for(j=0;j<1024;j++){
            double glat = (i - 256.0) * 90.0 / 256.0;
            double glon = -(j - 512.0) * 180.0 / 512.0;
            GalactictoRadec(glat, glon, &raa, &dec);
            double sindec = sin(dec);
            double cosdec = cos(dec);
            double sinlat = sin(lat);
            double coslat = cos(lat);
            radec_azel2(gstt - raa + lon, sindec, cosdec, sinlat, coslat, &azz, &el);
            fprintf(fp, "%.10e %.10e %.10e %.10e %.10e %.10e\n", glat, glon, raa, dec, azz, el);            
        }
    }
}