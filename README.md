# diffusion1
 diffusion1 - A package contining two pulse sequences for spin-echo
 diffusion
 measurements, "pge1.c" and "pge2.c". These sequence should be used
 with the Optional Diffusion Package where they supersede the standard
 sequences "pge.c" and "pgeramp.c". Both "pge1.c" and "pge2.c" measure
 diffusion with either a spin echo or a stimulated echo and allow for
 trapezoidal ramped gradients. "pge2.c" also provides a homospoil
 gradient during the stimulated-echo delay. Templates are provided
 for VNMR 6.1C and VnmrJ. Both sequences compile with Linux.

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

Submitter:      Dave Rice, Agilent

Date submitted: 2007-08-07 (VnmrJ 2.1B)
                2007-08-10 - Fixed bug, see below

File name:      diffusion1.tar.Z
Directory:      psglib

Description:    A package contining two pulse sequences for spin-echo
                diffusion measurements, "pge1.c" and "pge2.c". These sequence
                should be used with the Optional Diffusion Package where they
                supersede the standard sequences "pge.c" and "pgeramp.c". Both
                "pge1.c" and "pge2.c" measure diffusion with either a spin
                echo or a stimulated echo and allow for trapezoidal ramped
                gradients. "pge2.c" also provides a homospoil gradient during
                the stimulated-echo delay. Templates are provided for VNMR
                6.1C and VnmrJ. Both sequences compile with Linux.

                This package also supersedes the package "psglib/diffusion"
                which is only used with software up to VNMR 6.1C.

                Full use of these sequences requires the presence of the
                Optional Diffusion Package.

Related files:  Various files in execpars, maclib, manual, parlib, psglib,
                templates

Existing VnmrJ / VNMR files which are superseded or
otherwise affected by this submission:  None - supersedes "pge.c" and
                                           "pgeramp.c" of the optional
                                           diffusion package and updates
                                           "diffusion.tar.Z" for VnmrJ.
Hardware configuration limitations:     VNMRS (DirectDrive) ONLY
Known software version compatibility:      All versions up to VnmrJ 2.1B
Known OS version compatibility:         Linux - Solaris
Special instructions for installation:
    If you are downloading from the Internet, place the file "diffusion1.tar.Z"
    in "/vnmr/userlib/psglib" (you must be the primary VnmrJ user to do this).
    Log into the new user and type
        cd /vnmr/userlib
        ./extract psglib/diffusion1
        seqgen pge1.c pge2.c
    All files will load into the user's "~/vnmrsys" directory.

Reporting bugs:
---------------
Bug Reports are encouraged and should be reported through the standard bug
report at
        http://www.varianinc.com/products/nmr/apps/corner.html

**This software has not been tested on OpenVnmrJ. Use at your own risk.**

To install this user contribution:  
Download the repository from GitHub and checkout the tag of the contribution you want.
Typically tags end in the version (e.g. -v1.0)

     git clone https://github.com/OpenVnmrJ/psglib  
     cd psglib  
     git checkout diffusion1-v1.0


You may also download the archive directly from github at

    https://github.com/OpenVnmrJ/maclib/archive/diffusion1-v1.0.zip

Read diffusion1.README for installation instructions.

In most cases, move the contribution to /vnmr/userlib/psglib 
then use extract to install the contribution:  

    extract psglib/diffusion1