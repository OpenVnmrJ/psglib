# wpresat
 wpresat - a pulse sequence to do solvent suppression using obs.xmtr. The
 basic sequence is ..presat..p1..p1..at where all rf is waveform
 generator-based and all rf may be shaped

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

File name:      wpresat.c
Directory:      psglib
Description:    a pulse sequence to do solvent suppression using obs.xmtr.
                The basic sequence is ..presat..p1..p1..at where all rf
                is waveform generator-based and all rf may be shaped

Related files:  maclib/wpresat, manual/wpresat, parlib/wpresat.par,
                shapelib/presat

Existing VnmrJ / VNMR files which are superseded or
otherwise affected by this submission:  none
Hardware configuration limitations:     Requires linear amplifiers
                                        and waveform generator(1)
Known software version compatibility:   VnmrS 4.1
Known OS version compatibility:         n.a.
Special instructions for installation:

**This software has not been tested on OpenVnmrJ. Use at your own risk.**

To install this user contribution:  
Download the repository from GitHub and checkout the tag of the contribution you want.
Typically tags end in the version (e.g. -v1.0)

     git clone https://github.com/OpenVnmrJ/psglib  
     cd psglib  
     git checkout wpresat-v1.0


You may also download the archive directly from github at

    https://github.com/OpenVnmrJ/maclib/archive/wpresat-v1.0.zip

Read wpresat.README for installation instructions.

In most cases, move the contribution to /vnmr/userlib/psglib 
then use extract to install the contribution:  

    extract psglib/wpresat