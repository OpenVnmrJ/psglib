# gwsatnoe
 gwsatnoe - NOE difference experiment in which the irradition has
 multifrequency
 signal saturation using multifrequency SLP excitation. A shapelib
 file must exist that does a 90 degree soft pulse on all the desired
 lines simultaneously. This is followed by a gradient (homospoil)
 pulse to dephase the magnetization. This is repeated three times,
 increasing the phase of the pulse by 90 degrees each time. Finally,
 the whole four-step process is repeated multiple times.

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
Date submitted: 1995-02-09

File name:      gwsatnoe.c
Directory:      psglib
Description:    NOE difference experiment in which the irradition has
		multifrequency signal saturation using multifrequency SLP
		excitation. A shapelib file must exist that does a 90 degree
		soft pulse on all the desired lines simultaneously. This is
		followed by a gradient (homospoil) pulse to dephase the
		magnetization. This is repeated three times, increasing the
		phase of the pulse by 90 degrees each time. Finally, the whole
		four-step process is repeated multiple times.

		Create the soft pulse using the convolute macro (found in
		/vnmr/userlib). For example, a pulse "satpul" is created by:
                     convolute('gauss','satpul',10000,4300,2200,-124,-400)
		for saturating signals at resonance frequencies 4300,2200,
		-124, and -400 Hz relative to the xmtr position, using a
		gaussian pulse of 10msec duration (requires waveform gen-
		erator). On UNITYplus systems the pulse sequence may be
		modified to replace shaped_pulse with apshaped_pulse if no
		waveform generator is present in channel 1.

		If no pfg capability is present the sequence could be modified
		to use an ordinary homospoil.

		Sequence does noe difference by irradiating at dof and control
		on alternate scans and subtracting. nt should be a multiple of
		2 (large for better artifact suppression).

Related files:  parlib/gwsatnoe.par


Existing VnmrJ / VNMR files which are superseded or
otherwise affected by this submission:  none
Hardware configuration limitations:     pfg, wfg
Known software version compatibility:   VnmrS 4.3
Known OS version compatibility:         n.a.
Special instructions for installation:

**This software has not been tested on OpenVnmrJ. Use at your own risk.**

To install this user contribution:  
Download the repository from GitHub and checkout the tag of the contribution you want.
Typically tags end in the version (e.g. -v1.0)

     git clone https://github.com/OpenVnmrJ/psglib  
     cd psglib  
     git checkout gwsatnoe-v1.0


You may also download the archive directly from github at

    https://github.com/OpenVnmrJ/maclib/archive/gwsatnoe-v1.0.zip

Read gwsatnoe.README for installation instructions.

In most cases, move the contribution to /vnmr/userlib/psglib 
then use extract to install the contribution:  

    extract psglib/gwsatnoe