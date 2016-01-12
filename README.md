# nosyhmqc3d
 nosyhmqc3d - nosyhmqc3d with off-resonance presat (with shifted laminar
 pulse)
 option

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
Date submitted: 1992-10-28

File name:      nosyhmqc3d.c
Directory:      psglib
Description:    nosyhmqc3d with off-resonance presat (with shifted
                laminar pulse) option

Related files:  maclib/nosyhmqc3d, parlib/nosyhmqc3d13.par,
                parlib/nosyhmqc3d15.par, tablib/nosyhmqc3d, manual/nosyhmqc3d

Existing VnmrJ / VNMR files which are superseded or
otherwise affected by this submission:  earlier userlib version
Hardware configuration limitations:     PPM in channel 1
Known software version compatibility:   VnmrS 4.1A
Known OS version compatibility:         n.a.
Special instructions for installation:

**This software has not been tested on OpenVnmrJ. Use at your own risk.**

To install this user contribution:  
Download the repository from GitHub and checkout the tag of the contribution you want.
Typically tags end in the version (e.g. -v1.0)

     git clone https://github.com/OpenVnmrJ/psglib  
     cd psglib  
     git checkout nosyhmqc3d-v1.0


You may also download the archive directly from github at

    https://github.com/OpenVnmrJ/maclib/archive/nosyhmqc3d-v1.0.zip

Read nosyhmqc3d.README for installation instructions.

In most cases, move the contribution to /vnmr/userlib/psglib 
then use extract to install the contribution:  

    extract psglib/nosyhmqc3d