# goesytr
 goesytr - Double pulsed field gradient spinecho 1D ROE does transient 1D
 ROE
 experiment where very weak ROE's are visible because gradients dephase
 unwanted signals.  Selective pulses on desired multiplet can be easily
 created using setshape macro. Parameters are automatically determined
 and installed. The setshape macro can be run from ds_1 menu. No WFG
 is required when using VNMR 5.3 or later

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

Submitter:      George A. Gray, Varian
Date submitted: 1997-07-15
                1997-10-08 corrected use of "shape" parameter

File name:      goesytr
Directory:      psglib
Description:    Double pulsed field gradient spinecho 1D ROE does transient 1D
		ROE experiment where very weak ROE's are visible because
		gradients dephase unwanted signals.  Selective pulses on
		desired multiplet can be easily created using setshape macro.
		Parameters are automatically determined and installed.
		The setshape macro can be run from ds_1 menu.
                No WFG is required when using VNMR 5.3 or later

Related files:  menulib/ds_1,   maclib/goesytr,     maclib/setshape,
                manual/goesytr, parlib/goesytr.par, psglib/goesytr.c,

Existing VnmrJ / VNMR files which are superseded or
otherwise affected by this submission:  none
Hardware configuration limitations:     requires z-axis pfg
Known software version compatibility:   VNMR 5.3B
Known OS version compatibility:         n.a.
Special instructions for installation:
    If you are downloading from the Internet, store
    the file goesytr.tar.Z in /vnmr/userlib/psglib, then use
        cd /vnmr/userlib
        ./extract psglib/goesytr /vnmr
        seqgen goesytr
    for a "global" installation, for a local installation use
        cd /vnmr/userlib
        ./extract psglib/goesytr
        seqgen goesytr

**This software has not been tested on OpenVnmrJ. Use at your own risk.**

To install this user contribution:  
Download the repository from GitHub and checkout the tag of the contribution you want.
Typically tags end in the version (e.g. -v1.0)

     git clone https://github.com/OpenVnmrJ/psglib  
     cd psglib  
     git checkout goesytr-v1.0


You may also download the archive directly from github at

    https://github.com/OpenVnmrJ/maclib/archive/goesytr-v1.0.zip

Read goesytr.README for installation instructions.

In most cases, move the contribution to /vnmr/userlib/psglib 
then use extract to install the contribution:  

    extract psglib/goesytr