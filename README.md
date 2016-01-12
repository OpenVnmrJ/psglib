# hmqcnoesy3rf
 hmqcnoesy3rf - hmqcnoesy3rf uses third channel for N15 pulses and
 decoupling
 water suppression via  or presat replaces earlier userlib version
 hmqcnoesy3rf.c (fixed d2 corredtion to use pwx2 instead of pw)

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
Date submitted: 1993-01-13
		1994-01-31 (dec2pulse replaced with dec2rgpulse; TPPI
			    incrementation corrected; macro bug fixed - Krish)

File name:      hmqcnoesy3rf
Directory:      psglib
Description:    hmqcnoesy3rf uses third channel for N15 pulses and decoupling
                water suppression via  or presat
                replaces earlier userlib version hmqcnoesy3rf.c
                (fixed d2 corredtion to use pwx2 instead of pw)

Related files:  maclib/hmqcnoesy3rf, manual/hmqcnoesy3rf, parlib/hmqcnoesy3rf

Existing VnmrJ / VNMR files which are superseded or
otherwise affected by this submission:  none
Hardware configuration limitations:     third rf channel
Known software version compatibility:   VNMR 4.1, 4.3
Known OS version compatibility:         n.a.
Special instructions for installation:

**This software has not been tested on OpenVnmrJ. Use at your own risk.**

To install this user contribution:  
Download the repository from GitHub and checkout the tag of the contribution you want.
Typically tags end in the version (e.g. -v1.0)

     git clone https://github.com/OpenVnmrJ/psglib  
     cd psglib  
     git checkout hmqcnoesy3rf-v1.0


You may also download the archive directly from github at

    https://github.com/OpenVnmrJ/maclib/archive/hmqcnoesy3rf-v1.0.zip

Read hmqcnoesy3rf.README for installation instructions.

In most cases, move the contribution to /vnmr/userlib/psglib 
then use extract to install the contribution:  

    extract psglib/hmqcnoesy3rf