# The beamcorr algorithm as applied exactly to the Nature paper data.


```
# Known variables
integ = 6-18 hrs in seconds 
cmb = 2
Tcmb = 2.725

# Read Nivedit's beamfile (azelq.txt)
From this, we get:
frqspac = frq[1] - frq[0] of beam file. (= 2 MHz)
frqstart, freqstop as first and last frequencies in file between desired range.(here 40 and 100 MHz)

get azel as an array of the beam values, with shape (nfrq, naz, nel)
freqref = 75 MHz

nbeam = number of frequencies in azel array (31 here)

# READ THE SKY MODEL (408-all-noh.txt), which is a 512*1024 array covering the sky.
# We don't know the exact origin of this file. No corrections are applied in the code.

set lon = 116.5
set lat = -26.7  # both these values will be important.

set freq150 as the reference frequency, which is the one CLOSEST to freqref that's
actually in data.

get arrays of RA/DEC corresponding to the sky model, with shape (nra, ndec)


FOR each GHA in a grid from 6 to 18 with dGHA=0.5:
    >> For some reason, calculating the sun az/el at the gha (not sure anything is done with this)
    FOR dec in DEC array:
        get power-law model = (freq / 408 Mhz)**-2.5 (for freq in beam file)
        FOR ra in RA array:
            calculate az and el for this ra/dec
            get amp150 = \Omega_skymodel * beam(az, el, freq150)  # nearest-neighbors interpolation. remember freq150 is the closest frequency to freqref that's in the data.
            FOR freq IN beam freqs:
                amp = \Omega * beam(az, el, freq) # again, nearest-neighbors interpolation.
                if el is above horizon:
                    fg = (sky_model(ra, dec) - Tcmb) * power-law-model + Tcmb
                    bwfg += fg * amp        # fsum
                    bw += amp               # fsum1  
                    bwfg150 += fg * amp150  # fsum2

let Model = FourierModel(31 terms) # this is NOT smoothing

set beam_factor = (bwfg / bw) / (bwfg150(frq)  bwfg150(ref_freq))
fit the Model to the beam factor

```