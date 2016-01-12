# CIGAR
 CIGAR - Improved Accordion gradient HMBC. CIGAR2j3j does 2J/3J
 differentiation
 for protonated carbons.

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

Submitter:      Krish Krishnamurthy, Agilent
Date submitted: 2000-02-03
                2000-06-21 - added CIGAR2j3j

File name:      CIGAR
Directories:    psglib
Description:    Improved Accordion gradient HMBC. CIGAR2j3j does 2J/3J
                differentiation for protonated carbons.

Related files:  maclib, parlib, templates/dg, dialoglib, manual
                Also includes earlier CIGAR macros, sequence, etc.

Existing VnmrJ / VNMR files which are superseded or
otherwise affected by this submission:  earlier version of this contribution
Hardware configuration limitations:     Gradients
Known software version compatibility:   VNMR 6.1B
Known OS version compatibility:         n.a.
Special instructions for installation:
    If you are downloading from the Internet, store

**This software has not been tested on OpenVnmrJ. Use at your own risk.**

To install this user contribution:  
Download the repository from GitHub and checkout the tag of the contribution you want.
Typically tags end in the version (e.g. -v1.0)

     git clone https://github.com/OpenVnmrJ/psglib  
     cd psglib  
     git checkout CIGAR-v1.0


You may also download the archive directly from github at

    https://github.com/OpenVnmrJ/maclib/archive/CIGAR-v1.0.zip

Read CIGAR.README for installation instructions.

In most cases, move the contribution to /vnmr/userlib/psglib 
then use extract to install the contribution:  

    extract psglib/CIGAR