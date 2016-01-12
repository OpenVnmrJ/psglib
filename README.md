# 11echo
 11echo - water suppression by 11echo

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
Date submitted: 1991-12-10
                1994-01-31 (added rcvron and status(C);
                            hst reset to 2 ms in parameters - Krish)

File name:      11echo
Directory:      psglib
Description:    water suppression by 11echo

Related files:  maclib/11echo, manual/11echo, parlib/11echo


Existing VnmrJ / VNMR files which are superseded or
otherwise affected by this submission:  none
Hardware configuration limitations:     none
Known software version compatibility:   VNMR 4.1, 4.3
Known OS version compatibility:         n.a.
Special instructions for installation:

**This software has not been tested on OpenVnmrJ. Use at your own risk.**

To install this user contribution:  
Download the repository from GitHub and checkout the tag of the contribution you want.
Typically tags end in the version (e.g. -v1.0)

     git clone https://github.com/OpenVnmrJ/psglib  
     cd psglib  
     git checkout 11echo-v1.0


You may also download the archive directly from github at

    https://github.com/OpenVnmrJ/maclib/archive/11echo-v1.0.zip

Read 11echo.README for installation instructions.

In most cases, move the contribution to /vnmr/userlib/psglib 
then use extract to install the contribution:  

    extract psglib/11echo