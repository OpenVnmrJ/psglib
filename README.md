# gocsy
 gocsy - Double pulsed field gradient spinecho 1D tocsy Selective
 pulses on desired multiplet can be easily created using setshape
 macro. Parameters are automatically determined and installed. setshaped
 macro can be run from ds_1 menu. No WFG is required when using VNMR
 5.3 or later.

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

File name:      gocsy
Directory:      psglib
Description:    Double pulsed field gradient spinecho 1D tocsy Selective
		pulses on desired multiplet can be easily created using
		setshape macro. Parameters are automatically determined and
		installed. setshaped macro can be run from ds_1 menu.
                No WFG is required when using VNMR 5.3 or later.

Related files:  menulib/ds_1, maclib/gocsy,     maclib/setshape,
                manual/gocsy, parlib/gocsy.par, psglib/gocsy.c

Existing VnmrJ / VNMR files which are superseded or
otherwise affected by this submission:  none
Hardware configuration limitations:     requires z-axis pfg
Known software version compatibility:   VNMR 5.3B
Known OS version compatibility:         n.a.
Special instructions for installation:
    If you are downloading from the Internet, store
    the file gocsy.tar.Z in /vnmr/userlib/psglib, then use
        cd /vnmr/userlib
        ./extract psglib/gocsy /vnmr
        seqgen gocsy
    for a "global" installation, for a local installation use
        cd /vnmr/userlib
        ./extract psglib/gocsy
        seqgen gocsy

**This software has not been tested on OpenVnmrJ. Use at your own risk.**

To install this user contribution:  
Download the repository from GitHub and checkout the tag of the contribution you want.
Typically tags end in the version (e.g. -v1.0)

     git clone https://github.com/OpenVnmrJ/psglib  
     cd psglib  
     git checkout gocsy-v1.0


You may also download the archive directly from github at

    https://github.com/OpenVnmrJ/maclib/archive/gocsy-v1.0.zip

Read gocsy.README for installation instructions.

In most cases, move the contribution to /vnmr/userlib/psglib 
then use extract to install the contribution:  

    extract psglib/gocsy