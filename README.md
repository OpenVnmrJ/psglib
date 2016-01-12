# BIRDHMBC
 BIRDHMBC - BIRDHMQC (gradient HMBC with excellent suppression of residual
 1-bond
 correlations, superior to gHMBC, but 5 to 10% less sensitive than
 the latter), as published by Peter Bigler in Magn. Reson. Chem. 38,
 963-969 (2000).

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

Your name:              Bruce Adams
Address (optional):     9017 Mendenhall Ct., Suite D
                        Columbia, MD 21045
Phone (optional):       410-381-6024
FAX (optional):         410-381-7037
Date submitted:         2001-02-26

File name:              BIRDHMBC.c
Directory:              psglib
Description:            BIRDHMQC (gradient HMBC with excellent suppression of
                        residual 1-bond correlations, superior to gHMBC, but 5
                        to 10% less sensitive than the latter), as published by
                        Peter Bigler in Magn. Reson. Chem. 38, 963-969 (2000).

Related files:          dialoglib/BIRDHMBC, glide/exp/AuHexp/acquire.def,
                        maclib/BIRDHMBC, manual/BIRDHMBC, parlib/BIRDHMBC.par,
                        templates/dg/BIRDHMBC

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
     git checkout BIRDHMBC-v1.0


You may also download the archive directly from github at

    https://github.com/OpenVnmrJ/maclib/archive/BIRDHMBC-v1.0.zip

Read BIRDHMBC.README for installation instructions.

In most cases, move the contribution to /vnmr/userlib/psglib 
then use extract to install the contribution:  

    extract psglib/BIRDHMBC