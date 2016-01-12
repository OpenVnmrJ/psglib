# hmqctocsy3
 hmqctocsy3 - Replaces standard version in /vnmr/psglib with one which
 generates
 fewer acodes, and has reduced inter-increment delays

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

Your name:              Steve Patt
                        For support questions please contact
                                Rolf Kyburz (rolf.kyburz@agilent.com)
Date submitted: 1995-04-05

File name:      hmqctocsy3
Directory:      psglib
Description:    Replaces standard version in /vnmr/psglib with one which
                generates fewer acodes, and has reduced inter-increment
                delays

Related files:

Existing VnmrJ / VNMR files which are superseded or
otherwise affected by this submission:  psglib/hmqctocsy.c
Hardware configuration limitations:     Requires linear amplifiers
Known software version compatibility:   VnmrS 5.1
Known OS version compatibility:         n.a.
Special instructions for installation:

**This software has not been tested on OpenVnmrJ. Use at your own risk.**

To install this user contribution:  
Download the repository from GitHub and checkout the tag of the contribution you want.
Typically tags end in the version (e.g. -v1.0)

     git clone https://github.com/OpenVnmrJ/psglib  
     cd psglib  
     git checkout hmqctocsy3-v1.0


You may also download the archive directly from github at

    https://github.com/OpenVnmrJ/maclib/archive/hmqctocsy3-v1.0.zip

Read hmqctocsy3.README for installation instructions.

In most cases, move the contribution to /vnmr/userlib/psglib 
then use extract to install the contribution:  

    extract psglib/hmqctocsy3