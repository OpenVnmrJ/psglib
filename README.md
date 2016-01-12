# shphmqc
 shphmqc - hmqc using second channel for X pulses and decoupling uses shaped
 pulses for proton pulses

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
Date submitted: 1992-10-26
                1994-01-31 (macro & pulse sequence bugs fixed - Krish)

File name:      shphmqc
Directory:      psglib
Description:    hmqc using second channel for X pulses and decoupling
                uses shaped pulses for proton pulses

Related files:  maclib/shphmqc, manual/shphmqc, parlib/shphmqc

Existing VnmrJ / VNMR files which are superseded or
otherwise affected by this submission:  none
Hardware configuration limitations:     bb 2nd rf channel/wfg
Known software version compatibility:   VNMR 4.1, 4.3
Known OS version compatibility:         n.a.
Special instructions for installation:

**This software has not been tested on OpenVnmrJ. Use at your own risk.**

To install this user contribution:  
Download the repository from GitHub and checkout the tag of the contribution you want.
Typically tags end in the version (e.g. -v1.0)

     git clone https://github.com/OpenVnmrJ/psglib  
     cd psglib  
     git checkout shphmqc-v1.0


You may also download the archive directly from github at

    https://github.com/OpenVnmrJ/maclib/archive/shphmqc-v1.0.zip

Read shphmqc.README for installation instructions.

In most cases, move the contribution to /vnmr/userlib/psglib 
then use extract to install the contribution:  

    extract psglib/shphmqc