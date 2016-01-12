# wet
 wet - single or multiple line suppression using selective pulses followed
 by gradient homospoils. Can be used for water suppression and for
 multiple solvent suppression. Uses UNITYplus capability to control
 linear modulator to obtain fine-power control (can be modified to
 eliminate this for older systems). The pulse sequence may be used
 to replace first 90 degree pulse of other pulse sequences as form of
 solvent suppression. Selective pulses may be quite long (~20msec) and
 thus peaks very near the solvent can be observed with full intensity.

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
Date submitted: 1995-01-09
                1996-01-08 (modifies earlier version)

File name:      wet
Directory:      psglib
Description:    single or multiple line suppression using selective pulses
		followed by gradient homospoils. Can be used for water
		suppression and for multiple solvent suppression. Uses
		UNITYplus capability to control linear modulator to obtain
		fine-power control (can be modified to eliminate this for
		older systems). The pulse sequence may be used to replace
		first 90 degree pulse of other pulse sequences as form of
		solvent suppression. Selective pulses may be quite long
		(~20msec) and thus peaks very near the solvent can be observed
		with full intensity.

		Organic solvents having multiple lines, or mixtures of
		protonated solvents can be suppressed by having the shaped
		pulse saturate all simultaneously (see gwsatnoe.REAME for how
		to do this). C13 satellites of these organic solvents will
		possibly be still very strong. These may be eliminated by
		decoupling C13 during the shaped pulse using "dm"

                Sequence written by Steve Smallcombe, Varian

Related files:  parlib/wet.par tablib/wet


Existing VnmrJ / VNMR files which are superseded or
otherwise affected by this submission:  Earlier version in userlib which
					didn't reset fine power before final
					pulse
Hardware configuration limitations:     pfg
Known software version compatibility:   VNMR 5.1
Known OS version compatibility:         n.a.
Special instructions for installation:

**This software has not been tested on OpenVnmrJ. Use at your own risk.**

To install this user contribution:  
Download the repository from GitHub and checkout the tag of the contribution you want.
Typically tags end in the version (e.g. -v1.0)

     git clone https://github.com/OpenVnmrJ/psglib  
     cd psglib  
     git checkout wet-v1.0


You may also download the archive directly from github at

    https://github.com/OpenVnmrJ/maclib/archive/wet-v1.0.zip

Read wet.README for installation instructions.

In most cases, move the contribution to /vnmr/userlib/psglib 
then use extract to install the contribution:  

    extract psglib/wet