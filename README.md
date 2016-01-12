# wetseq
 wetseq - This packages contains a series of 1D and 2D pulse sequences
 incorporating WET solvent suppression for liquids NMR as described
 in J.Magn.Reson. A 117, 295 (1995). WET uses shaped pulses and pulsed
 field gradients for suppression of one or more solvents.

 Copyright 2016 University of Oregon

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

                                SUBMISSION FORM

Your name:      Steve Smallcombe, Varian / Paul Keifer, Agilent
Date submitted: 1996-03-14
                1996-07-01 (added lcpar2d, lcpsgset, lcset2d macros)

File name:      wet
Directory:      psglib
Description:    This packages contains a series of 1D and 2D pulse sequences
		incorporating WET solvent suppression for liquids NMR as
		described in J.Magn.Reson. A 117, 295 (1995). WET uses shaped
		pulses and pulsed field gradients for suppression of one or
		more solvents.

Related files:  psglib/wet1d.c       psglib/wetdqcosy.c   psglib/wetgcosy.c
                psglib/wetghmqcps.c  psglib/wetghsqc.c    psglib/wetgmqcosyps.c
                psglib/wetnoesy.c    psglib/wetpwxcal.c   psglib/wetrelayh.c
                psglib/wettntocsy.c  shapelib/seduce1.RF  shapelib/wet.RF
                and related files.

Existing VnmrJ / VNMR files which are superseded or
otherwise affected by this submission:  none
Hardware configuration limitations:     PFG, wfg
        UNITYplus, UNITY INOVA; UNITY systems require PSG fix, see below.
Known software version compatibility:   VNMR 5.1
Special instructions for installation:
    If you are downloading from the Internet, store
    the file wetseq.tar.Z in /vnmr/userlib/psglib, then use
        cd /vnmr/userlib
        ./extract psglib/wetseq /vnmr
        seqgen ~/vnmrsys/psglib/wet*.c    <seqgen /vnmr/psglib/wet*.c>
    for a "global" installation, for a local installation use
        cd /vnmr/userlib
        ./extract psglib/wetseq
        seqgen ~/vnmrsys/psglib/wet*.c    <seqgen /vnmr/psglib/wet*.c>
    The WET pulse sequence element takes advantage of the scaling of a shaped
    pulse using fine power adjustment (pwrf function). In VNMR 5.1 this
    functionality is disabled for UNITY systems, to avoid a conflict with
    systems using a fine attenuator. It is possible to make tpwf scaling work
    on UNITY systems, using the following procedure:
        setuserpsg
        cd ~/vnmrsys/psg
        cp /vnmr/psg/wg.c .

**This software has not been tested on OpenVnmrJ. Use at your own risk.**

To install this user contribution:  
Download the repository from GitHub and checkout the tag of the contribution you want.
Typically tags end in the version (e.g. -v1.0)

     git clone https://github.com/OpenVnmrJ/psglib  
     cd psglib  
     git checkout wetseq-v1.0


You may also download the archive directly from github at

    https://github.com/OpenVnmrJ/maclib/archive/wetseq-v1.0.zip

Read wetseq.README for installation instructions.

In most cases, move the contribution to /vnmr/userlib/psglib 
then use extract to install the contribution:  

    extract psglib/wetseq