#Userlib/psglib
This repository contains the userlib/psglib contributions for VnmrJ as
packaged in VnmrJ 4.2

##NOTES & CAUTIONS
* Many of these contributions were incorporated into the core VnmrJ
software years ago, so may be redundant with core capabilities
* Most of the contributions are obsolete.
* Though many contributions work, there is often a better way of doing
things in the more modern software
* These contributions are not guaranteed to work
* These tools may work on some systems but not others. Many will only
work with certain versions of VnmrJ, VNMR, RHEL, or Solaris.
* Use these tools at your own risk. Some of them, such as pulse sequences,
could potentially damage hardware, and neither Agilent nor UO is
responsible for such damage.
* Though this is the User Library, users' contributions are mixed
with those from Varian and Agilent staff. On Spinsights, similar
Agilent-provided materials like Chempack are available in Toolkit, and
shared materials from users/customers are here in User Library. Agilent
staff still do contribute to the currently active User Library, but,
like all other contributions, those are personal materials that they
have developed and find useful, and they are not officially supported by
Agilent or guaranteed to work.
* The last update to this file was performed on February 1, 2013, and it
will not be updated again.
* This initially should only contain contributions from Agilent and/or
Varian but your contributions are welcome

##Downloading
You may download this repository from GitHub at:
https://github.com/OpenVnmrJ/psglib.git

##Updating and adding
- Fork on GitHub
- Do not add your contribution to the master branch
- If updating contribution, checkout the tag of the contribution, update
and commit on the contribution branch
- If adding a new contribution, checkout a new branch with the name of
the contribution, push the new branch to your repository, add and commit
on the new branch
- Tag your update of change with the name of the contribution, followed
by a version, for example, mymacro-v1.1
- Push your branch to your fork; remember to push the tags too
- Make a pull request to the OpenVnmrJ repository

##Contributions

Below is a list of each contribution. To access a contributions, check
out its tag

##11echo
>water suppression by 11echo

---
To install the contribution, checkout the tag 11echo-v1.0:

    git checkout 11echo-v1.0

then read 11echo.README

Usually use extract to install the contribution:

    extract psglib/11echo

##adequate_AD
>Pulse sequence for 1-1 Adequate, uses adiabatic and "av180" 13C 180
degree pulses. References: - B.Rief, M.Koeck, R.Kerssebaum, H.Kang,
W.Fenical, and C.Griesinger, J. Magn. Reson. A 118, 282 (1996). -
M.Koeck, B.Rief, W.Fenical, and C.Griesinger, Tetrahedron Lett. 37(3),
363 (1996).

---
To install the contribution, checkout the tag adequate_AD-v1.0:

    git checkout adequate_AD-v1.0

then read adequate_AD.README

Usually use extract to install the contribution:

    extract psglib/adequate_AD

##bashdtoxy
>BAnd Selective Homonuclear Decoupled TOXY (Ref.: V.V. Krishnamurthy,
Magn. Reson. Chem., in press)

---
To install the contribution, checkout the tag bashdtoxy-v1.0:

    git checkout bashdtoxy-v1.0

then read bashdtoxy.README

Usually use extract to install the contribution:

    extract psglib/bashdtoxy

##BIRDHMBC
>BIRDHMQC (gradient HMBC with excellent suppression of residual 1-bond
correlations, superior to gHMBC, but 5 to 10% less sensitive than
the latter), as published by Peter Bigler in Magn. Reson. Chem. 38,
963-969 (2000).

---
To install the contribution, checkout the tag BIRDHMBC-v1.0:

    git checkout BIRDHMBC-v1.0

then read BIRDHMBC.README

Usually use extract to install the contribution:

    extract psglib/BIRDHMBC

##caconh
>caconh uses third channel for N15 pulses and decoupling; uses second
channel for 13C pulses; water suppression via gradients and purge
pulses; uses phase-cycling for coherence selection. Sequence submitted
and written by Gordon Rule, U.Virginia and moving to Carnegie-Mellon
Univ. Spring 1996)

---
To install the contribution, checkout the tag caconh-v1.0:

    git checkout caconh-v1.0

then read caconh.README

Usually use extract to install the contribution:

    extract psglib/caconh

##cconh
>This contribution includes the following pulse sequences (based on
code written by Lewis Kay, U.Toronto), including auxiliary macros for
the creation of the necessary SLP pulse shapes: - c_co_nh , known as
C(CO)NH in the literature; - two gradient versions of C(CO)NH, with
and without seduce decoupling in place of a 180 C13 pulse, which
appear to give equivalent results; - the equivalent two gradient
versions of CBCA(CO)NH.

---
To install the contribution, checkout the tag cconh-v1.0:

    git checkout cconh-v1.0

then read cconh.README

Usually use extract to install the contribution:

    extract psglib/cconh

##cctocsy
>13C Tocsy followed by inept to detected 1H

---
To install the contribution, checkout the tag cctocsy-v1.0:

    git checkout cctocsy-v1.0

then read cctocsy.README

Usually use extract to install the contribution:

    extract psglib/cctocsy

##CIGAR
>Improved Accordion gradient HMBC. CIGAR2j3j does 2J/3J differentiation
for protonated carbons.

---
To install the contribution, checkout the tag CIGAR-v1.0:

    git checkout CIGAR-v1.0

then read CIGAR.README

Usually use extract to install the contribution:

    extract psglib/CIGAR

##coconoesy
>Concurent cosy/noesy  experiment (J.Magn.Reson. 56, 343(1984)). Two
versions are supplied, coconoesy.c & coconoesyT.c, both of which
support a random mix time. coconoesyT.c is a tablib-based random delay
version, whereas coconoesy.c has the random delay generated from the
older "randomdelay.c" routine. If you desire random delays the "T"
version definately works better.  Either version allows one to use
a fixed delay.

---
To install the contribution, checkout the tag coconoesy-v1.0:

    git checkout coconoesy-v1.0

then read coconoesy.README

Usually use extract to install the contribution:

    extract psglib/coconoesy

##cpmgt2
>carr-purcell meiboom-gill T2 experiment

---
To install the contribution, checkout the tag cpmgt2-v1.0:

    git checkout cpmgt2-v1.0

then read cpmgt2.README

Usually use extract to install the contribution:

    extract psglib/cpmgt2

##cycledof
>noe-difference via multiplet cycling using decoupler for saturation

---
To install the contribution, checkout the tag cycledof-v1.0:

    git checkout cycledof-v1.0

then read cycledof.README

Usually use extract to install the contribution:

    extract psglib/cycledof

##cyclenoePPM
>noe-difference via multiplet cycling using the programmable pulse
modulator (PPM) for selective excitation and control pulses. May be
used for obs xmtr irradiation and attenuation of selective excitation
and control pulses. eg.  if shape=hard102.RF such that have 1/10
amplitude of shape=hard.RF then a 20 db attentuation of a standard
hard pulse can be provided. Related files:    maclib/cyclenoePPM,
manual/cyclenoePPM, parlib/cyclenoePPM, shapelib/hard102.RF

---
To install the contribution, checkout the tag cyclenoePPM-v1.0:

    git checkout cyclenoePPM-v1.0

then read cyclenoePPM.README

Usually use extract to install the contribution:

    extract psglib/cyclenoePPM

##cyclenoeSLP
>noe-difference via multiplet cycling using obs channel for saturation
and pulsing and using SLP pulses

---
To install the contribution, checkout the tag cyclenoeSLP-v1.0:

    git checkout cyclenoeSLP-v1.0

then read cyclenoeSLP.README

Usually use extract to install the contribution:

    extract psglib/cyclenoeSLP

##cyclenoe
>NOE-difference via multiplet cycling using the observe channel for
saturation and pulsing (this sequence is part of the standard "psglib"
in VnmrJ 2.2C & later versions).

---
To install the contribution, checkout the tag cyclenoe-v1.0:

    git checkout cyclenoe-v1.0

then read cyclenoe.README

Usually use extract to install the contribution:

    extract psglib/cyclenoe

##dantesat
>presat using DANTE off-resonance suppression

---
To install the contribution, checkout the tag dantesat-v1.0:

    git checkout dantesat-v1.0

then read dantesat.README

Usually use extract to install the contribution:

    extract psglib/dantesat

##dante
>DANTE off-resonance excitation: (....pulse)n

---
To install the contribution, checkout the tag dante-v1.0:

    git checkout dante-v1.0

then read dante.README

Usually use extract to install the contribution:

    extract psglib/dante

##deptghmqc
>depthmqc with gradient selection

---
To install the contribution, checkout the tag deptghmqc-v1.0:

    git checkout deptghmqc-v1.0

then read deptghmqc.README

Usually use extract to install the contribution:

    extract psglib/deptghmqc

##deptgl
>optimized dept

---
To install the contribution, checkout the tag deptgl-v1.0:

    git checkout deptgl-v1.0

then read deptgl.README

Usually use extract to install the contribution:

    extract psglib/deptgl

##depthmqc
>depthmqc experiment

---
To install the contribution, checkout the tag depthmqc-v1.0:

    git checkout depthmqc-v1.0

then read depthmqc.README

Usually use extract to install the contribution:

    extract psglib/depthmqc

##depttocsy
>depthmqc variant of hmqctocsy

---
To install the contribution, checkout the tag depttocsy-v1.0:

    git checkout depttocsy-v1.0

then read depttocsy.README

Usually use extract to install the contribution:

    extract psglib/depttocsy

##diffusion1
>A package contining two pulse sequences for spin-echo diffusion
measurements, "pge1.c" and "pge2.c". These sequence should be used
with the Optional Diffusion Package where they supersede the standard
sequences "pge.c" and "pgeramp.c". Both "pge1.c" and "pge2.c" measure
diffusion with either a spin echo or a stimulated echo and allow for
trapezoidal ramped gradients. "pge2.c" also provides a homospoil
gradient during the stimulated-echo delay. Templates are provided
for VNMR 6.1C and VnmrJ. Both sequences compile with Linux.

---
To install the contribution, checkout the tag diffusion1-v1.0:

    git checkout diffusion1-v1.0

then read diffusion1.README

Usually use extract to install the contribution:

    extract psglib/diffusion1

##dipshft
>sequence, parlib and manual entry for simple dipshift (SLF) (corrected
9/15/92 to eliminate "gate" function)

---
To install the contribution, checkout the tag dipshft-v1.0:

    git checkout dipshft-v1.0

then read dipshft.README

Usually use extract to install the contribution:

    extract psglib/dipshft

##ecosy
>"Exclusive" COSY (E.COSY)

---
To install the contribution, checkout the tag ecosy-v1.0:

    git checkout ecosy-v1.0

then read ecosy.README

Usually use extract to install the contribution:

    extract psglib/ecosy

##exside
>EXcitation Sculptured Indirect Detection Experiment with J-scaling
option (V.V. Krishnamurthy, JMR, A121, 33, 1996)

---
To install the contribution, checkout the tag exside-v1.0:

    git checkout exside-v1.0

then read exside.README

Usually use extract to install the contribution:

    extract psglib/exside

##fhoesy
>F19-detected 2D proton fluorine NOE experiment

---
To install the contribution, checkout the tag fhoesy-v1.0:

    git checkout fhoesy-v1.0

then read fhoesy.README

Usually use extract to install the contribution:

    extract psglib/fhoesy

##flock
>flock  experiment

---
To install the contribution, checkout the tag flock-v1.0:

    git checkout flock-v1.0

then read flock.README

Usually use extract to install the contribution:

    extract psglib/flock

##ghmqcps
>hmqc using gradients for coherence pathway seletion with
phase-sensitive detection in F1 (Does not do gradient HMBC!)

---
To install the contribution, checkout the tag ghmqcps-v1.0:

    git checkout ghmqcps-v1.0

then read ghmqcps.README

Usually use extract to install the contribution:

    extract psglib/ghmqcps

##ghsqcse
>hsqc using gradients for coherence pathway seletion. Uses
Rance,Cavanaugh sensitivity enhancement. Uses gradient pulses for
eliminating imperfections in 180's. Has "flipback" soft pulse for
reducing H2O saturation. See Kay, Keifer and Saarinen, JACS, 114,
10663(1992)

---
To install the contribution, checkout the tag ghsqcse-v1.0:

    git checkout ghsqcse-v1.0

then read ghsqcse.README

Usually use extract to install the contribution:

    extract psglib/ghsqcse

##gmqcosyps
>cosy using gradients for coherence pathway seletion with
phase-sensitive detection in F1

---
To install the contribution, checkout the tag gmqcosyps-v1.0:

    git checkout gmqcosyps-v1.0

then read gmqcosyps.README

Usually use extract to install the contribution:

    extract psglib/gmqcosyps

##gocsy
>Double pulsed field gradient spinecho 1D tocsy Selective
pulses on desired multiplet can be easily created using setshape
macro. Parameters are automatically determined and installed. setshaped
macro can be run from ds_1 menu. No WFG is required when using VNMR
5.3 or later.

---
To install the contribution, checkout the tag gocsy-v1.0:

    git checkout gocsy-v1.0

then read gocsy.README

Usually use extract to install the contribution:

    extract psglib/gocsy

##goesy
>Double pulsed field gradient spinecho NOE does transient 1D NOE
experiment where very weak NOE's are visible because gradients dephase
unwanted signals.  Selective pulses on desired multiplet can be easily
created using setshape macro. Parameters are automatically determined
and installed. The setshape macro can be run from ds_1 menu. No WFG
is required when using VNMR 5.3 or later

---
To install the contribution, checkout the tag goesy-v1.0:

    git checkout goesy-v1.0

then read goesy.README

Usually use extract to install the contribution:

    extract psglib/goesy

##goesytr
>Double pulsed field gradient spinecho 1D ROE does transient 1D ROE
experiment where very weak ROE's are visible because gradients dephase
unwanted signals.  Selective pulses on desired multiplet can be easily
created using setshape macro. Parameters are automatically determined
and installed. The setshape macro can be run from ds_1 menu. No WFG
is required when using VNMR 5.3 or later

---
To install the contribution, checkout the tag goesytr-v1.0:

    git checkout goesytr-v1.0

then read goesytr.README

Usually use extract to install the contribution:

    extract psglib/goesytr

##gsh2pul
>s2pul with shaped gradients

---
To install the contribution, checkout the tag gsh2pul-v1.0:

    git checkout gsh2pul-v1.0

then read gsh2pul.README

Usually use extract to install the contribution:

    extract psglib/gsh2pul

##gtocsy
>tocsy using optional non-water excitation

---
To install the contribution, checkout the tag gtocsy-v1.0:

    git checkout gtocsy-v1.0

then read gtocsy.README

Usually use extract to install the contribution:

    extract psglib/gtocsy

##gwsatnoe
>NOE difference experiment in which the irradition has multifrequency
signal saturation using multifrequency SLP excitation. A shapelib
file must exist that does a 90 degree soft pulse on all the desired
lines simultaneously. This is followed by a gradient (homospoil)
pulse to dephase the magnetization. This is repeated three times,
increasing the phase of the pulse by 90 degrees each time. Finally,
the whole four-step process is repeated multiple times.

---
To install the contribution, checkout the tag gwsatnoe-v1.0:

    git checkout gwsatnoe-v1.0

then read gwsatnoe.README

Usually use extract to install the contribution:

    extract psglib/gwsatnoe

##hacaconh
>caconh uses third channel for N15 pulses and decoupling; uses second
channel for 13C pulses; water suppression via gradients, and purge
pulse; uses phase-cycling for coherence selection. (sequence submitted
and written by Gordon Rule, U.Virginia and moving to Carnegie-Mellon
Univ. Spring 1996)

---
To install the contribution, checkout the tag hacaconh-v1.0:

    git checkout hacaconh-v1.0

then read hacaconh.README

Usually use extract to install the contribution:

    extract psglib/hacaconh

##hahnse
>Hahn spin-echo as per E. L. Hahn, Phys. Rev. 80, 580-594, 1950. Also
used for solvent/H2O suppression as per Rabenstein et al., Analytical
Chem. 60, 13 80A, 1988.

---
To install the contribution, checkout the tag hahnse-v1.0:

    git checkout hahnse-v1.0

then read hahnse.README

Usually use extract to install the contribution:

    extract psglib/hahnse

##hcacon_ct
>hcacon using third channel for N15 pulses and decoupling; uses
second channel for 13C decoupling. Constant time version (adapted
from code by S.Smallcombe); corrected previous version (needed +6db
on cshape3 power)

---
To install the contribution, checkout the tag hcacon_ct-v1.0:

    git checkout hcacon_ct-v1.0

then read hcacon_ct.README

Usually use extract to install the contribution:

    extract psglib/hcacon_ct

##hcannh_HC
>hcannh_HC uses third channel for N15 pulses and decoupling uses second
channel for 13C decoupling 4D pulse sequence with F1=alpha protons,
F2=alpha carbons

---
To install the contribution, checkout the tag hcannh_HC-v1.0:

    git checkout hcannh_HC-v1.0

then read hcannh_HC.README

Usually use extract to install the contribution:

    extract psglib/hcannh_HC

##hcannh_NC
>hcannh_NC uses third channel for N15 pulses and decoupling uses second
channel for 13C decoupling 4D pulse sequence with F1=NH protons,
F2=alpha carbons

---
To install the contribution, checkout the tag hcannh_NC-v1.0:

    git checkout hcannh_NC-v1.0

then read hcannh_NC.README

Usually use extract to install the contribution:

    extract psglib/hcannh_NC

##hcchcosy
>a pulse sequence to do relate protons via CC cosy in 13C-labeled
molecules. Modified from code by Wittekamp and Mueller.

---
To install the contribution, checkout the tag hcchcosy-v1.0:

    git checkout hcchcosy-v1.0

then read hcchcosy.README

Usually use extract to install the contribution:

    extract psglib/hcchcosy

##hcchtocsy2
>Replaces standard version in /vnmr/psglib with one which generates
fewer acodes, and has reduced inter-increment delays

---
To install the contribution, checkout the tag hcchtocsy2-v1.0:

    git checkout hcchtocsy2-v1.0

then read hcchtocsy2.README

Usually use extract to install the contribution:

    extract psglib/hcchtocsy2

##hcchtocsy
>A pulse sequence to do relate protons via CC tocsy in 13C-labeled
molecules. Modified from code by Wittekamp and Mueller.

---
To install the contribution, checkout the tag hcchtocsy-v1.0:

    git checkout hcchtocsy-v1.0

then read hcchtocsy.README

Usually use extract to install the contribution:

    extract psglib/hcchtocsy

##hcosyinv
>heteronuclear cosy with proton observe

---
To install the contribution, checkout the tag hcosyinv-v1.0:

    git checkout hcosyinv-v1.0

then read hcosyinv.README

Usually use extract to install the contribution:

    extract psglib/hcosyinv

##hetcorcp
>Heteronuclear correlation in the solid-state

---
To install the contribution, checkout the tag hetcorcp-v1.0:

    git checkout hetcorcp-v1.0

then read hetcorcp.README

Usually use extract to install the contribution:

    extract psglib/hetcorcp

##hetcorps
>phase-sensitive and absolute value hetcor with option to do
multiplicity selection; improved from earlier User Library version

---
To install the contribution, checkout the tag hetcorps-v1.0:

    git checkout hetcorps-v1.0

then read hetcorps.README

Usually use extract to install the contribution:

    extract psglib/hetcorps

##heteroHH
>heteronuclear Hartmann-Hahn X-H correlation uses DIPSI mixing and
channel 2 or 3

---
To install the contribution, checkout the tag heteroHH-v1.0:

    git checkout heteroHH-v1.0

then read heteroHH.README

Usually use extract to install the contribution:

    extract psglib/heteroHH

##hetsecho_3rf
>heteronuclear spin-echo difference experiment using second decoupler

---
To install the contribution, checkout the tag hetsecho_3rf-v1.0:

    git checkout hetsecho_3rf-v1.0

then read hetsecho_3rf.README

Usually use extract to install the contribution:

    extract psglib/hetsecho_3rf

##hetsecho
>heteronuclear spin-echo difference experiment

---
To install the contribution, checkout the tag hetsecho-v1.0:

    git checkout hetsecho-v1.0

then read hetsecho.README

Usually use extract to install the contribution:

    extract psglib/hetsecho

##hmqc11_3rf
>hmqc using third channel for X pulses and decoupling uses 11 pulses
for 90 and 180

---
To install the contribution, checkout the tag hmqc11_3rf-v1.0:

    git checkout hmqc11_3rf-v1.0

then read hmqc11_3rf.README

Usually use extract to install the contribution:

    extract psglib/hmqc11_3rf

##hmqc3rf_C13dec
>hmqc using third channel for X pulses and decoupling with the option
of decoupling C13 during t1 (for use in doubly-enriched molecules)

---
To install the contribution, checkout the tag hmqc3rf_C13dec-v1.0:

    git checkout hmqc3rf_C13dec-v1.0

then read hmqc3rf_C13dec.README

Usually use extract to install the contribution:

    extract psglib/hmqc3rf_C13dec

##hmqc3rf
>hmqc using third channel for X pulses and decoupling (corrected
version of previous hmqc3rf.c)

---
To install the contribution, checkout the tag hmqc3rf-v1.0:

    git checkout hmqc3rf-v1.0

then read hmqc3rf.README

Usually use extract to install the contribution:

    extract psglib/hmqc3rf

##hmqc_COdec
>hmqc using second channel for X pulses and decoupling with the option
of decoupling 13CO and 15N during t1 (for use in doubly-enriched
molecules)

---
To install the contribution, checkout the tag hmqc_COdec-v1.0:

    git checkout hmqc_COdec-v1.0

then read hmqc_COdec.README

Usually use extract to install the contribution:

    extract psglib/hmqc_COdec

##hmqcjr_3rf
>hmqc using third channel for X pulses and decoupling uses jr(jump
and return) pulses for 90 and 180 (corrected pulse sequence from 10
Dec 1991 version) (aka 11echo version of hmqc - as per Bax)

---
To install the contribution, checkout the tag hmqcjr_3rf-v1.0:

    git checkout hmqcjr_3rf-v1.0

then read hmqcjr_3rf.README

Usually use extract to install the contribution:

    extract psglib/hmqcjr_3rf

##hmqcnoesy3rf
>hmqcnoesy3rf uses third channel for N15 pulses and decoupling
water suppression via  or presat replaces earlier userlib version
hmqcnoesy3rf.c (fixed d2 corredtion to use pwx2 instead of pw)

---
To install the contribution, checkout the tag hmqcnoesy3rf-v1.0:

    git checkout hmqcnoesy3rf-v1.0

then read hmqcnoesy3rf.README

Usually use extract to install the contribution:

    extract psglib/hmqcnoesy3rf

##hmqcnoesy
>hmqc followed by noesy (F1=X, F2=1H)

---
To install the contribution, checkout the tag hmqcnoesy-v1.0:

    git checkout hmqcnoesy-v1.0

then read hmqcnoesy.README

Usually use extract to install the contribution:

    extract psglib/hmqcnoesy

##hmqc
>hmqc from VnmrS 4.1 modified to add "sspul" and to correct d2 for
precession during pwx, plus minor changes in parameters

---
To install the contribution, checkout the tag hmqc-v1.0:

    git checkout hmqc-v1.0

then read hmqc.README

Usually use extract to install the contribution:

    extract psglib/hmqc

##hmqctocsy2
>hmqctocsy: HYPERCOMPL.FAD with editing of direct responses

---
To install the contribution, checkout the tag hmqctocsy2-v1.0:

    git checkout hmqctocsy2-v1.0

then read hmqctocsy2.README

Usually use extract to install the contribution:

    extract psglib/hmqctocsy2

##hmqctocsy3
>Replaces standard version in /vnmr/psglib with one which generates
fewer acodes, and has reduced inter-increment delays

---
To install the contribution, checkout the tag hmqctocsy3-v1.0:

    git checkout hmqctocsy3-v1.0

then read hmqctocsy3.README

Usually use extract to install the contribution:

    extract psglib/hmqctocsy3

##hmqctocsyS
>hmqctocsy with F1 selective options. FAD + direct response editing.

---
To install the contribution, checkout the tag hmqctocsyS-v1.0:

    git checkout hmqctocsyS-v1.0

then read hmqctocsyS.README

Usually use extract to install the contribution:

    extract psglib/hmqctocsyS

##hmqctocsy
>hmqc followed by tocsy (F1=X, F2=1H)

---
To install the contribution, checkout the tag hmqctocsy-v1.0:

    git checkout hmqctocsy-v1.0

then read hmqctocsy.README

Usually use extract to install the contribution:

    extract psglib/hmqctocsy

##hmqctoxy3d3rf
>hmqc followed by tocsy (F1=X, F2=F3=1H) (corrected previous version
-removed gate function)

---
To install the contribution, checkout the tag hmqctoxy3d3rf-v1.0:

    git checkout hmqctoxy3d3rf-v1.0

then read hmqctoxy3d3rf.README

Usually use extract to install the contribution:

    extract psglib/hmqctoxy3d3rf

##hmqctoxy3d
>HMQC-TOCSY-3D experiment. Improved axial suppression.

---
To install the contribution, checkout the tag hmqctoxy3d-v1.0:

    git checkout hmqctoxy3d-v1.0

then read hmqctoxy3d.README

Usually use extract to install the contribution:

    extract psglib/hmqctoxy3d

##HMSC
>Simultaneously detected heteronuclear shift correlation through
multiple and single bonds (the latter is an absolute value, coupled
HMQC). Ref.: R. Burger, C. Schorn, & P.  Bigler, J. Magn. Reson., 148,
88-94 (2001). For direct correlation result use wft2d(1,0,1,0,0,1,0,1)
for long-range correlation use wft2d(1,0,-1,0,0,1,0,-1)

---
To install the contribution, checkout the tag HMSC-v1.0:

    git checkout HMSC-v1.0

then read HMSC.README

Usually use extract to install the contribution:

    extract psglib/HMSC

##hncanodec
>H-->N-->C(alpha) correlation experiment Farmer, Venters, Spicer,
Wittekind and Mueller, J.Magn. Reson. (in press)

---
To install the contribution, checkout the tag hncanodec-v1.0:

    git checkout hncanodec-v1.0

then read hncanodec.README

Usually use extract to install the contribution:

    extract psglib/hncanodec

##hnca
>H-->N-->C(alpha) correlation experiment Farmer, Venters, Spicer,
Wittekind and Mueller, J.Magn. Reson. (in press) (changed hnca macro)

---
To install the contribution, checkout the tag hnca-v1.0:

    git checkout hnca-v1.0

then read hnca.README

Usually use extract to install the contribution:

    extract psglib/hnca

##hncoca_alt
>H-->N-->C=O-->C(alpha) correlation experiment

---
To install the contribution, checkout the tag hncoca_alt-v1.0:

    git checkout hncoca_alt-v1.0

then read hncoca_alt.README

Usually use extract to install the contribution:

    extract psglib/hncoca_alt

##hncocanodec
>H-->N-->C=O-->C(alpha) correlation experiment

---
To install the contribution, checkout the tag hncocanodec-v1.0:

    git checkout hncocanodec-v1.0

then read hncocanodec.README

Usually use extract to install the contribution:

    extract psglib/hncocanodec

##hncoca_simple
>H-->N-->C=O-->C(alpha) correlation experiment THIS IS A PULSE SEQUENCE
THAT HAS BEEN SIMPLIFIED FOR USE BY "DPS". IT IS NOT AS ACCURATE AS
HNCOCA AND SHOULD NOT BE USED TO ACCUMULATE DATA. IT WILL, HOWEVER,
BE USEFUL FOR DISPLAY OF THE PULSE SEQUENCE

---
To install the contribution, checkout the tag hncoca_simple-v1.0:

    git checkout hncoca_simple-v1.0

then read hncoca_simple.README

Usually use extract to install the contribution:

    extract psglib/hncoca_simple

##hncoca
>hncoca uses third channel for N15 pulses and decoupling uses second
channel for 13C pulses water suppression via trim pulse or presat
replaces earlier userlib version hncoca_simple.c (this is based on
the original published version by bax and ikura, but uses 3 channels
instead of 4)

---
To install the contribution, checkout the tag hncoca-v1.0:

    git checkout hncoca-v1.0

then read hncoca.README

Usually use extract to install the contribution:

    extract psglib/hncoca

##hnconodec
>H-->N-->C=O correlation experiment

---
To install the contribution, checkout the tag hnconodec-v1.0:

    git checkout hnconodec-v1.0

then read hnconodec.README

Usually use extract to install the contribution:

    extract psglib/hnconodec

##hnco
>H-->N-->C=O correlation experiment

---
To install the contribution, checkout the tag hnco-v1.0:

    git checkout hnco-v1.0

then read hnco.README

Usually use extract to install the contribution:

    extract psglib/hnco

##HNHSS3D
>hmqc-noesy-hmqc-3D sequence using SS pulses for final hmqc step to
avoid presaturaton. Uses 3rd channel for X.

---
To install the contribution, checkout the tag HNHSS3D-v1.0:

    git checkout HNHSS3D-v1.0

then read HNHSS3D.README

Usually use extract to install the contribution:

    extract psglib/HNHSS3D

##hs90
>hs90 sequence required for CRAMPS/Multipulse installation.

---
To install the contribution, checkout the tag hs90-v1.0:

    git checkout hs90-v1.0

then read hs90.README

Usually use extract to install the contribution:

    extract psglib/hs90

##hsqc3rf_C13dec
>hsqc using third channel for X pulses and decoupling uses second
channel for optional 13C decoupling in di-labelled molecules. uses
status-controlled C=O decoupling in t1 with C(alpha) decoupling using
an SLP 180 degree pulse in t1.

---
To install the contribution, checkout the tag hsqc3rf_C13dec-v1.0:

    git checkout hsqc3rf_C13dec-v1.0

then read hsqc3rf_C13dec.README

Usually use extract to install the contribution:

    extract psglib/hsqc3rf_C13dec

##hsqcT1
>hsqcT1 uses third channel for N15 pulses and decoupling uses second
channel for 13C decoupling(under dm control) water suppression via
trim pulse or presat proton decoupling during t1

---
To install the contribution, checkout the tag hsqcT1-v1.0:

    git checkout hsqcT1-v1.0

then read hsqcT1.README

Usually use extract to install the contribution:

    extract psglib/hsqcT1

##hsqc
>Heteronuclear Single-Quantum Correlation

---
To install the contribution, checkout the tag hsqc-v1.0:

    git checkout hsqc-v1.0

then read hsqc.README

Usually use extract to install the contribution:

    extract psglib/hsqc

##hsqctoxySE
>hsqc-tocsy-3D experiment with Rance method sensitivity enhancement
in BOTH t1 and t2 dimensions

---
To install the contribution, checkout the tag hsqctoxySE-v1.0:

    git checkout hsqctoxySE-v1.0

then read hsqctoxySE.README

Usually use extract to install the contribution:

    extract psglib/hsqctoxySE

##jsdqcosy_table
>j-scaled double-quantum-filtered cosy using phase tables

---
To install the contribution, checkout the tag jsdqcosy_table-v1.0:

    git checkout jsdqcosy_table-v1.0

then read jsdqcosy_table.README

Usually use extract to install the contribution:

    extract psglib/jsdqcosy_table

##jsdqcosy
>j-scaled double-quantum-filtered cosy

---
To install the contribution, checkout the tag jsdqcosy-v1.0:

    git checkout jsdqcosy-v1.0

then read jsdqcosy.README

Usually use extract to install the contribution:

    extract psglib/jsdqcosy

##jumpret
>water suppression by jump-and-return

---
To install the contribution, checkout the tag jumpret-v1.0:

    git checkout jumpret-v1.0

then read jumpret.README

Usually use extract to install the contribution:

    extract psglib/jumpret

##lrhetcor
>long-range heteronuclear chemical shift correlation

---
To install the contribution, checkout the tag lrhetcor-v1.0:

    git checkout lrhetcor-v1.0

then read lrhetcor.README

Usually use extract to install the contribution:

    extract psglib/lrhetcor

##noesy11
>noesy using 11-echo for last pulse

---
To install the contribution, checkout the tag noesy11-v1.0:

    git checkout noesy11-v1.0

then read noesy11.README

Usually use extract to install the contribution:

    extract psglib/noesy11

##noesyhmqc3rf
>noesy followed by hmqc filter (F1=X, F2=1H) X pulses and decoupling
uses third channel

---
To install the contribution, checkout the tag noesyhmqc3rf-v1.0:

    git checkout noesyhmqc3rf-v1.0

then read noesyhmqc3rf.README

Usually use extract to install the contribution:

    extract psglib/noesyhmqc3rf

##noesyjr45
>noesy with jump-return read pulse phases for phase=1 and phase=2 are
-45 and 45 degrees to balance residual water signal evenly in both
data tables (removed rcvroff/rcvron from previous version)

---
To install the contribution, checkout the tag noesyjr45-v1.0:

    git checkout noesyjr45-v1.0

then read noesyjr45.README

Usually use extract to install the contribution:

    extract psglib/noesyjr45

##noesyjr
>noesy with jumpr-return read pulse

---
To install the contribution, checkout the tag noesyjr-v1.0:

    git checkout noesyjr-v1.0

then read noesyjr.README

Usually use extract to install the contribution:

    extract psglib/noesyjr

##noesyskbx
>noesy using a sklenar-bax monitoring "pulse"

---
To install the contribution, checkout the tag noesyskbx-v1.0:

    git checkout noesyskbx-v1.0

then read noesyskbx.README

Usually use extract to install the contribution:

    extract psglib/noesyskbx

##noesysl
>noesy with trim -free precession - trim uses no presat. trim
pulses dephase water. free precession period spreads signals as in
jump-return(same amplitude dependence)

---
To install the contribution, checkout the tag noesysl-v1.0:

    git checkout noesysl-v1.0

then read noesysl.README

Usually use extract to install the contribution:

    extract psglib/noesysl

##noesywg_r
>noesywg with suppression of radiation damping in t1 and mix

---
To install the contribution, checkout the tag noesywg_r-v1.0:

    git checkout noesywg_r-v1.0

then read noesywg_r.README

Usually use extract to install the contribution:

    extract psglib/noesywg_r

##nosyhmqc3d3rf
>noesy followed by hmqc filter (F1=F3=1H, F2=X)

---
To install the contribution, checkout the tag nosyhmqc3d3rf-v1.0:

    git checkout nosyhmqc3d3rf-v1.0

then read nosyhmqc3d3rf.README

Usually use extract to install the contribution:

    extract psglib/nosyhmqc3d3rf

##nosyhmqc3d
>nosyhmqc3d with off-resonance presat (with shifted laminar pulse)
option

---
To install the contribution, checkout the tag nosyhmqc3d-v1.0:

    git checkout nosyhmqc3d-v1.0

then read nosyhmqc3d.README

Usually use extract to install the contribution:

    extract psglib/nosyhmqc3d

##nosytoxy3d
>noesy - tocsy - 3D

---
To install the contribution, checkout the tag nosytoxy3d-v1.0:

    git checkout nosytoxy3d-v1.0

then read nosytoxy3d.README

Usually use extract to install the contribution:

    extract psglib/nosytoxy3d

##nosyxfil
>noesy with w2 x-filter with off-resonance presat option

---
To install the contribution, checkout the tag nosyxfil-v1.0:

    git checkout nosyxfil-v1.0

then read nosyxfil.README

Usually use extract to install the contribution:

    extract psglib/nosyxfil

##presati
>a pulse sequence to do solvent suppression using obs.xmtr. in which
pulse sequence elements are active in IPA (in ACQI)

---
To install the contribution, checkout the tag presati-v1.0:

    git checkout presati-v1.0

then read presati.README

Usually use extract to install the contribution:

    extract psglib/presati

##presat
>a pulse sequence to do solvent suppression using obs.xmtr.

---
To install the contribution, checkout the tag presat-v1.0:

    git checkout presat-v1.0

then read presat.README

Usually use extract to install the contribution:

    extract psglib/presat

##refhsqc1
>a pulse sequence to do H--X correlations using a refocussed
double-inverse INEPT. It is written to allow X operation on either
channel of a 3channel spectrometer It permits proton decoupling during
t1. Differs from refhsqc by reduced # of pulses Author: Sandy Farmer,
Varian

---
To install the contribution, checkout the tag refhsqc1-v1.0:

    git checkout refhsqc1-v1.0

then read refhsqc1.README

Usually use extract to install the contribution:

    extract psglib/refhsqc1

##refhsqc
>Heteronuclear Single-Quantum Correlation

---
To install the contribution, checkout the tag refhsqc-v1.0:

    git checkout refhsqc-v1.0

then read refhsqc.README

Usually use extract to install the contribution:

    extract psglib/refhsqc

##roesy1d
>1D roesy using waveform generator for selective pulse

---
To install the contribution, checkout the tag roesy1d-v1.0:

    git checkout roesy1d-v1.0

then read roesy1d.README

Usually use extract to install the contribution:

    extract psglib/roesy1d

##sanc1
>Offset alternating sequence as per Lee et. al. "Multipole Theory of
Composite Pulses II: Time Reversal and Offset Alternating Sequences
With Experimental Confirmation", Journal of Magnetic Resonance 92,
455, 1991

---
To install the contribution, checkout the tag sanc1-v1.0:

    git checkout sanc1-v1.0

then read sanc1.README

Usually use extract to install the contribution:

    extract psglib/sanc1

##selexcit
>DPFGSE band-selection 1D

---
To install the contribution, checkout the tag selexcit-v1.0:

    git checkout selexcit-v1.0

then read selexcit.README

Usually use extract to install the contribution:

    extract psglib/selexcit

##selhmqc3rf
>selhmqc using third channel for X pulses and decoupling

---
To install the contribution, checkout the tag selhmqc3rf-v1.0:

    git checkout selhmqc3rf-v1.0

then read selhmqc3rf.README

Usually use extract to install the contribution:

    extract psglib/selhmqc3rf

##selhmqc
>selective hmqc

---
To install the contribution, checkout the tag selhmqc-v1.0:

    git checkout selhmqc-v1.0

then read selhmqc.README

Usually use extract to install the contribution:

    extract psglib/selhmqc

##shphmqc_3rf
>hmqc using third channel for X pulses and decoupling uses shaped
pulses for proton pulses

---
To install the contribution, checkout the tag shphmqc_3rf-v1.0:

    git checkout shphmqc_3rf-v1.0

then read shphmqc_3rf.README

Usually use extract to install the contribution:

    extract psglib/shphmqc_3rf

##shphmqc
>hmqc using second channel for X pulses and decoupling uses shaped
pulses for proton pulses

---
To install the contribution, checkout the tag shphmqc-v1.0:

    git checkout shphmqc-v1.0

then read shphmqc.README

Usually use extract to install the contribution:

    extract psglib/shphmqc

##simba
>hmqc/hmbc with selective options added

---
To install the contribution, checkout the tag simba-v1.0:

    git checkout simba-v1.0

then read simba.README

Usually use extract to install the contribution:

    extract psglib/simba

##sklbax
>a pulse sequence for solvent suppression using obs.xmtr. only The
basic sequence is ..d1..softH2O90.hard90&trim..at where the soft
pulse is waveform generator-based and may be shaped (wsklbax is the
same but requires waveform gen.)

---
To install the contribution, checkout the tag sklbax-v1.0:

    git checkout sklbax-v1.0

then read sklbax.README

Usually use extract to install the contribution:

    extract psglib/sklbax

##slpsat
>A variety of standard pulse sequences which use the SLP technique
(shifted laminar pulses) to produce off-resonance presaturation.

---
To install the contribution, checkout the tag slpsat-v1.0:

    git checkout slpsat-v1.0

then read slpsat.README

Usually use extract to install the contribution:

    extract psglib/slpsat

##sshoesy1d
>1D steady-state heteronuclear NOE experiment. Ref:  C. Hofstetter &
T.C. Pochapsky, Magn. Reson. Chem. 38, 90-94 (2000).

---
To install the contribution, checkout the tag sshoesy1d-v1.0:

    git checkout sshoesy1d-v1.0

then read sshoesy1d.README

Usually use extract to install the contribution:

    extract psglib/sshoesy1d

##SSnoesy
>noesy pulse sequence with SS read pulse

---
To install the contribution, checkout the tag SSnoesy-v1.0:

    git checkout SSnoesy-v1.0

then read SSnoesy.README

Usually use extract to install the contribution:

    extract psglib/SSnoesy

##superw
>A solvent suppression pulse sequence which is self-refocusing for
excited spins.

---
To install the contribution, checkout the tag superw-v1.0:

    git checkout superw-v1.0

then read superw.README

Usually use extract to install the contribution:

    extract psglib/superw

##swarp3d
>3D spinwarp imaging sequence

---
To install the contribution, checkout the tag swarp3d-v1.0:

    git checkout swarp3d-v1.0

then read swarp3d.README

Usually use extract to install the contribution:

    extract psglib/swarp3d

##t1cp
>T1 measurement in solids via cross polarization-method of D.A. Torchia,
JMR, v30, 613, (1978)

---
To install the contribution, checkout the tag t1cp-v1.0:

    git checkout t1cp-v1.0

then read t1cp.README

Usually use extract to install the contribution:

    extract psglib/t1cp

##tncosyps
>a pulse sequence to do cosyps with water suppression and pulses coming
from only the  obs.xmtr.

---
To install the contribution, checkout the tag tncosyps-v1.0:

    git checkout tncosyps-v1.0

then read tncosyps.README

Usually use extract to install the contribution:

    extract psglib/tncosyps

##tndqcosy
>a pulse sequence to do dqcosy with water suppression and pulses coming
from only the  obs.xmtr.  (changed macro)

---
To install the contribution, checkout the tag tndqcosy-v1.0:

    git checkout tndqcosy-v1.0

then read tndqcosy.README

Usually use extract to install the contribution:

    extract psglib/tndqcosy

##tnmqcosy
>multiple-quantum filtered cosy. Observe channel used for presaturation

---
To install the contribution, checkout the tag tnmqcosy-v1.0:

    git checkout tnmqcosy-v1.0

then read tnmqcosy.README

Usually use extract to install the contribution:

    extract psglib/tnmqcosy

##tnnoesy
>a pulse sequence to do noesy with water suppression and pulses coming
from only the  obs.xmtr.

---
To install the contribution, checkout the tag tnnoesy-v1.0:

    git checkout tnnoesy-v1.0

then read tnnoesy.README

Usually use extract to install the contribution:

    extract psglib/tnnoesy

##tnroesy
>a pulse sequence to do roesy with water suppression and pulses coming
from only the  obs.xmtr. corrected previous version: satmode control
in relax. delay

---
To install the contribution, checkout the tag tnroesy-v1.0:

    git checkout tnroesy-v1.0

then read tnroesy.README

Usually use extract to install the contribution:

    extract psglib/tnroesy

##tntocsy
>a pulse sequence to do tocsy with water suppression and pulses coming
from only the  obs.xmtr. (corrected version without "gate" function)

---
To install the contribution, checkout the tag tntocsy-v1.0:

    git checkout tntocsy-v1.0

then read tntocsy.README

Usually use extract to install the contribution:

    extract psglib/tntocsy

##tocsy1d1
>a pulse sequence to do 1D tocsy with water suppression and pulses
coming from only the  obs.xmtr.. The selective inversion pulse is
generated via a apshaped_pulse pulse sequence element that ramps the
linear modulator

---
To install the contribution, checkout the tag tocsy1d1-v1.0:

    git checkout tocsy1d1-v1.0

then read tocsy1d1.README

Usually use extract to install the contribution:

    extract psglib/tocsy1d1

##tocsy1d
>a pulse sequence to do 1D tocsy with water suppression and pulses
coming from only the  obs.xmtr.. The selective inversion pulse is
generated via a attn_shape_pulse pulse sequence element that ramps the
programmable attenuator through half-gauss, hermite or gaussian shapes.

---
To install the contribution, checkout the tag tocsy1d-v1.0:

    git checkout tocsy1d-v1.0

then read tocsy1d.README

Usually use extract to install the contribution:

    extract psglib/tocsy1d

##tqfb
>Triple Quantum Filter for measuring biexponential transverse relaxation
times of Quadrupolar nuclei (Na23 Li7);As per Chun-Wa Chung and
S. Wimperis, J. Magn. Reson. 88, 440, 1990.

---
To install the contribution, checkout the tag tqfb-v1.0:

    git checkout tqfb-v1.0

then read tqfb.README

Usually use extract to install the contribution:

    extract psglib/tqfb

##trnsySS3d
>3D sequence to do TOCSY-NOESY-3D, TOCSY-ROESY-3D, ROESY-NOESY-3D,
with SS read pulse

---
To install the contribution, checkout the tag trnsySS3d-v1.0:

    git checkout trnsySS3d-v1.0

then read trnsySS3d.README

Usually use extract to install the contribution:

    extract psglib/trnsySS3d

##troesy1
>transverse cross-relaxation experiment in rotating frame (ref: Shaka,
JACS 114 3157 (1992))

---
To install the contribution, checkout the tag troesy1-v1.0:

    git checkout troesy1-v1.0

then read troesy1.README

Usually use extract to install the contribution:

    extract psglib/troesy1

##w1dztocsy_slp
>a pulse sequence to do 1dtocsy with water suppression and pulses coming
from only the obs.xmtr. The selective inversion pulse is generated
via a shaped_pulse pulse sequence element that uses the waveform
generator. Z-filtering is an option. The frequency remains at "tof"
during the whole experiment. Pulses and inversion pulses must be
"shifted laminar pulses" if off-resonance. corrected to eliminate
"gate" functions

---
To install the contribution, checkout the tag w1dztocsy_slp-v1.0:

    git checkout w1dztocsy_slp-v1.0

then read w1dztocsy_slp.README

Usually use extract to install the contribution:

    extract psglib/w1dztocsy_slp

##w1dztocsy
>a pulse sequence to do 1dtocsy with water suppression and pulses
coming from only the obs.xmtr. The selective inversion pulse is
generated via a shaped_pulse pulse sequence element that uses the
waveform generator. Z-filtering is an option (corrected to eliminate
"gate" functions)

---
To install the contribution, checkout the tag w1dztocsy-v1.0:

    git checkout w1dztocsy-v1.0

then read w1dztocsy.README

Usually use extract to install the contribution:

    extract psglib/w1dztocsy

##watergnoe
>NOESY with water suppression via WATERGATE gradient echo

---
To install the contribution, checkout the tag watergnoe-v1.0:

    git checkout watergnoe-v1.0

then read watergnoe.README

Usually use extract to install the contribution:

    extract psglib/watergnoe

##waterg
>WATERGATE 1D with water suppression via gradient echo

---
To install the contribution, checkout the tag waterg-v1.0:

    git checkout waterg-v1.0

then read waterg.README

Usually use extract to install the contribution:

    extract psglib/waterg

##wetseq
>This packages contains a series of 1D and 2D pulse sequences
incorporating WET solvent suppression for liquids NMR as described
in J.Magn.Reson. A 117, 295 (1995). WET uses shaped pulses and pulsed
field gradients for suppression of one or more solvents.

---
To install the contribution, checkout the tag wetseq-v1.0:

    git checkout wetseq-v1.0

then read wetseq.README

Usually use extract to install the contribution:

    extract psglib/wetseq

##wet
>single or multiple line suppression using selective pulses followed
by gradient homospoils. Can be used for water suppression and for
multiple solvent suppression. Uses UNITYplus capability to control
linear modulator to obtain fine-power control (can be modified to
eliminate this for older systems). The pulse sequence may be used
to replace first 90 degree pulse of other pulse sequences as form of
solvent suppression. Selective pulses may be quite long (~20msec) and
thus peaks very near the solvent can be observed with full intensity.

---
To install the contribution, checkout the tag wet-v1.0:

    git checkout wet-v1.0

then read wet.README

Usually use extract to install the contribution:

    extract psglib/wet

##wlexch
>2D wideline exchange experiment (Spiess method)

---
To install the contribution, checkout the tag wlexch-v1.0:

    git checkout wlexch-v1.0

then read wlexch.README

Usually use extract to install the contribution:

    extract psglib/wlexch

##wnoesy
>a pulse sequence to do noesy using obs.xmtr. for all rf. The basic
sequence is ..presat..p1.d2.p2.mix.p3..at where all rf is waveform
generator-based and all rf may be shaped Selective noesy is done by
making p1 a selective pulse and reducing F1 (the p1 pulse may be a
"shifted laminar pulse" to effect a 90 degree pulse at a specific
position away from the transmitter- the shapelib entry listed below
is one for a 90 degree tophat pulse done 1080 Hz away from tof).

---
To install the contribution, checkout the tag wnoesy-v1.0:

    git checkout wnoesy-v1.0

then read wnoesy.README

Usually use extract to install the contribution:

    extract psglib/wnoesy

##wpresat
>a pulse sequence to do solvent suppression using obs.xmtr. The
basic sequence is ..presat..p1..p1..at where all rf is waveform
generator-based and all rf may be shaped

---
To install the contribution, checkout the tag wpresat-v1.0:

    git checkout wpresat-v1.0

then read wpresat.README

Usually use extract to install the contribution:

    extract psglib/wpresat

##wsklbax
>a pulse sequence for solvent suppression using obs.xmtr. only The
basic sequence is ..d1..softH2O90.hard90&trim..at where the soft
pulse is waveform generator-based and may be shaped

---
To install the contribution, checkout the tag wsklbax-v1.0:

    git checkout wsklbax-v1.0

then read wsklbax.README

Usually use extract to install the contribution:

    extract psglib/wsklbax

##xdec2pul
>A pulse sequence to perform F19 decouple - 1H and F19-1H NOE difference
experiments on single broadband UNITY spectrometers (fixed frequency
decoupler).

---
To install the contribution, checkout the tag xdec2pul-v1.0:

    git checkout xdec2pul-v1.0

then read xdec2pul.README

Usually use extract to install the contribution:

    extract psglib/xdec2pul

