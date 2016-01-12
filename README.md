# refhsqc1
 refhsqc1 - a pulse sequence to do H--X correlations using a refocussed
 double-inverse INEPT. It is written to allow X operation on either
 channel of a 3channel spectrometer It permits proton decoupling during
 t1. Differs from refhsqc by reduced # of pulses Author: Sandy Farmer,
 Varian

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
Date submitted: 1991-10-07
                1993-01-08 (removed rlpower definition - r.k.)
                1994-01-31 (added macro, parameters - Krish)

File name:      refhsqc1
Directory:      psglib
Description:    a pulse sequence to do H--X correlations
                using a refocussed double-inverse INEPT. It is written to
                allow X operation on either channel of a 3channel spectrometer
                It permits proton decoupling during t1.
                Differs from refhsqc by reduced # of pulses
		Author: Sandy Farmer, Varian

Related files:  maclib/refhscq1, manual/refhscq1, parlib/refhscq1

Existing VnmrJ / VNMR files which are superseded or
otherwise affected by this submission:  none
Hardware configuration limitations:     Requires linear amplifiers
Known software version compatibility:   VNMR 4.1, 4.3
Known OS version compatibility:         n.a.
Special instructions for installation:

**This software has not been tested on OpenVnmrJ. Use at your own risk.**

To install this user contribution:  
Download the repository from GitHub and checkout the tag of the contribution you want.
Typically tags end in the version (e.g. -v1.0)

     git clone https://github.com/OpenVnmrJ/psglib  
     cd psglib  
     git checkout refhsqc1-v1.0


You may also download the archive directly from github at

    https://github.com/OpenVnmrJ/maclib/archive/refhsqc1-v1.0.zip

Read refhsqc1.README for installation instructions.

In most cases, move the contribution to /vnmr/userlib/psglib 
then use extract to install the contribution:  

    extract psglib/refhsqc1