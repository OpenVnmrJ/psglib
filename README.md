# HMSC
 HMSC - Simultaneously detected heteronuclear shift correlation through
 multiple and single bonds (the latter is an absolute value, coupled
 HMQC). Ref.: R. Burger, C. Schorn, & P.  Bigler, J. Magn. Reson., 148,
 88-94 (2001). For direct correlation result use wft2d(1,0,1,0,0,1,0,1)
 for long-range correlation use wft2d(1,0,-1,0,0,1,0,-1)

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

Your name:      Bruce Adams, Varian
Date submitted: 2001-02-26

File name:      HMSC.c
Directory:      psglib
Description:    Simultaneously detected heteronuclear shift correlation
		through multiple and single bonds (the latter is an absolute
		value, coupled HMQC). Ref.: R. Burger, C. Schorn, & P.  Bigler,
		J. Magn. Reson., 148, 88-94 (2001).
                For direct correlation result use
                        wft2d(1,0,1,0,0,1,0,1)
                for long-range correlation use
                        wft2d(1,0,-1,0,0,1,0,-1)

Related files:  dialoglib/HMSC, glide/exp/AuHexp/acquire.def,
                maclib/HMSC, manual/HMSC, parlib/HMSC.par, templates/dg/HMSC

Existing VnmrJ / VNMR files which are superseded or
otherwise affected by this submission:  glide/exp/AuHexp/acquire.def
Hardware configuration limitations:     PFG and linear amplifiers
Known software version compatibility:   VNMR 6.1B
Known OS version compatibility:         Solaris 2.6
Special instructions for installation:

**This software has not been tested on OpenVnmrJ. Use at your own risk.**

To install this user contribution:  
Download the repository from GitHub and checkout the tag of the contribution you want.
Typically tags end in the version (e.g. -v1.0)

     git clone https://github.com/OpenVnmrJ/psglib  
     cd psglib  
     git checkout HMSC-v1.0


You may also download the archive directly from github at

    https://github.com/OpenVnmrJ/maclib/archive/HMSC-v1.0.zip

Read HMSC.README for installation instructions.

In most cases, move the contribution to /vnmr/userlib/psglib 
then use extract to install the contribution:  

    extract psglib/HMSC