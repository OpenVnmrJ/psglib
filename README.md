# cyclenoe
 cyclenoe - NOE-difference via multiplet cycling using the observe channel
 for
 saturation and pulsing (this sequence is part of the standard "psglib"
 in VnmrJ 2.2C & later versions).

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

Submitter:      George A. Gray, Agilent

Date submitted: 1991-12-10
                1994-01-31 "abort()" changed to "abort(1)" - Krish
                1994-02-02 Implemented Steve McKenna's proposed fix for short
                           saturation times - r.k.
                2011-01-18 Updated for current hard- & software - Paul Keifer

File name:      cyclenoe
Directory:      psglib
Description:    NOE-difference via multiplet cycling using the observe channel
                for saturation and pulsing (this sequence is part of the
                standard "psglib" in VnmrJ 2.2C & later versions).

Related files:  maclib/cyclenoe, manual/cyclenoe, parlib/cyclenoe

Existing VnmrJ / VNMR files which are superseded or
otherwise affected by this submission:  none
Hardware configuration limitations:     none with current instruments
Known software version compatibility:   VnmrJ 2.1 (no VnmrJ panels provided)
Known OS version compatibility:         n.a.
Special instructions for installation:

**This software has not been tested on OpenVnmrJ. Use at your own risk.**

To install this user contribution:  
Download the repository from GitHub and checkout the tag of the contribution you want.
Typically tags end in the version (e.g. -v1.0)

     git clone https://github.com/OpenVnmrJ/psglib  
     cd psglib  
     git checkout cyclenoe-v1.0


You may also download the archive directly from github at

    https://github.com/OpenVnmrJ/maclib/archive/cyclenoe-v1.0.zip

Read cyclenoe.README for installation instructions.

In most cases, move the contribution to /vnmr/userlib/psglib 
then use extract to install the contribution:  

    extract psglib/cyclenoe