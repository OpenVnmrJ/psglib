/* cyclenoeSLP. Difference NOE experiment
      Uses the obs xmtr for irradiation
      and frequency-shifted rf (for single-frequency irradiation)
      using shape "satshape" can be shaped and frequency-shifted,
      or shaped using "offset" to shift, depending on values of satfrq and tof

      frequency-cycled irradiation only done via offset PSG element
  Parameters:

    pw - 90 degree excitation pulse (at power tpwr)
    intsub - 'y': internal subtraction of data acquired by on-resonance
                  and off-resonance selective excitation
             'n': data acquired by on-resonance and off-resonance selec-
                  tive excitation are stored separately
    satfrq  - frequency of selective saturation (on resonance);
              (make this an array if intsub = 'n', off-res dof if cycle='n')
    control - off resonance selective saturation frequency
              (an inactive parameter if intsub = 'n')
    cycle   - cycle='y' does on-resonance saturation using frequency cycling
              around the frequency "satfrq" given by "spacing"and"pattern" 
              cycle='n' does off-resonance saturation at control
    spacing - spacing(in Hz) of multiplet
    pattern - pattern type ( 1 for singlet, 2 for doublet, etc.)    
    tau     - period spent on a single irradiation point during cycling 
    satpwr  - power of selective irradiation
    sattime - total length of irradiation at frequency satfrq.
  satshape  - name of shapelib file for saturation pulse
  conshape  - name of shapelib file for control pulse
    mix - mixing time
    hs - 'ynyn'
    sspul   - sspul='y' does trim(x)-trim(y)  before d1
    nt - intsub = 'n':  multiple of 16
         intsub = 'y':  multiple of 32


    NOTE:  This pulse sequence requires that the observe channel be
           equipped with direct synthesis RF and a linear amplifier.
           satpwr ranges 0-63 control obs attenuator
       Contact-  G.Gray (palo alto)    revision- from cycledof.c*/


#include <standard.h>

pulsesequence()
{
/*DEFINE LOCAL VARIABLES */
  double mix,satfrq,control,satpwr,sattime,
         tau,spacing;
  int pattern,times,jj;
  char conshape[MAXSTR],satshape[MAXSTR],intsub[MAXSTR],cycle[MAXSTR],sspul[MAXSTR];


/* LOAD AND INITIALIZE VARIABLES */
  tau = getval("tau"); getstr("satshape",satshape);
  getstr("conshape",conshape);
  getstr("intsub",intsub); getstr("cycle",cycle); getstr("sspul",sspul); 
  satfrq = getval("satfrq"); control = getval("control");
  sattime = getval("sattime"); spacing = getval("spacing");
  satpwr = getval("satpwr"); pattern = (int)getval("pattern");
  mix = getval("mix");
  if (pattern == 0) pattern = 1; if (tau == 0.0) tau = 0.1;
  times = (int)(sattime/(pattern*tau));
    initval(satpwr,v6);
    initval(0.0,v10);
  initval(0.0,v5);
  initval(tpwr,v4);

/* CHECK CONDITIONS */
 if (!newtrans)
   {
    fprintf(stdout,"REQUIRED:  direct syn. RF and linear amplifiers.\n");
    abort(1);
   }


/* CALCULATE PHASES */
 if (intsub[0] == 'y') hlv(ct,v1); else assign(ct,v1);
 assign(v1,oph);
 if (intsub[0] == 'y')
   {
    mod2(ct,v14);  /* trigger for the alteration of the saturation freq */
    ifzero(v14); add(oph,two,oph); endif(v14);
   }


/* BEGIN ACTUAL PULSE SEQUENCE CODE */
 status(A);
    if (sspul[0] == 'y')
    {
   pulse(1000*pw,zero); pulse(1000*pw,one);;
    }
    hsdelay(d1);
    power(v6,TODEV);
 status(B);

/* selective pulse or decoupler saturation */
    /* no cycling or interleaved subtraction (make control an array)*/
    if ((intsub[0] == 'n') && (cycle[0] == 'n'))
     {
      if (fabs(tof - satfrq) >=0.1) offset(satfrq,TODEV);
      shaped_pulse(satshape,sattime,zero,rof1,rof2);
      if (fabs(tof - satfrq) >=0.1) offset(tof,TODEV);
     }
    /* interleaved subtraction without cycling */
    if ((intsub[0] == 'y') && (cycle[0] == 'n'))
      {
       ifzero(v14);
          shaped_pulse(conshape,sattime,zero,rof1,rof2);
       elsenz(v14);
          if (fabs(tof - satfrq) >=0.1) offset(satfrq,TODEV);
          shaped_pulse(satshape,sattime,zero,rof1,rof2);
          if (fabs(tof - satfrq) >=0.1) offset(tof,TODEV);
       endif(v14);
      }
    /* no interleaved subtraction but cycling is used (make cycle array)*/
    if ((cycle[0] == 'y') && (intsub[0] == 'n'))
    {
     for (jj = 0; jj < times; jj++)
      {
       double startfrq; int i;
       startfrq = satfrq - (pattern/2)*spacing;
       if ((pattern %2) == 0) startfrq = startfrq + spacing/2.0;
       for (i = 0; i < pattern; i++)
        {
         offset(startfrq,TODEV);
         rgpulse(tau,zero,rof1,rof2);
         startfrq = startfrq + spacing;
        }
      }  
    }
    /* interleaved subtraction with cycling (no array needed for one
       value of satfrq. Link array satfrq with pattern and spacing for
       multiple noe difference spectra within one experiment. For example
       set array = '(satfrq,pattern,spacing)') */

    if ((cycle[0] == 'y') && (intsub[0] == 'y'))
     {
      ifzero(v14);
       for (jj = 0; jj < times; jj++)
        {
         double startfrq; int i;
         startfrq = control - (pattern/2)*spacing;
         if ((pattern %2) == 0) startfrq = startfrq + spacing/2.0;
         for (i = 0; i < pattern; i++)
          {
           offset(startfrq,TODEV);
           rgpulse(tau,zero,rof1,rof2);
           startfrq = startfrq + spacing;
          }
         }
      elsenz(v14);
       for (jj = 0; jj < times; jj++)
        {
         double startfrq; int i;
         startfrq = satfrq - (pattern/2)*spacing;
         if ((pattern %2) == 0) startfrq = startfrq + spacing/2.0;
         for (i = 0; i < pattern; i++)
          {
           offset(startfrq,TODEV);
           rgpulse(tau,zero,rof1,rof2);
           startfrq = startfrq + spacing;
          }
        }  
      endif(v14);
     }
    /* restore power levels as controlled by tpwr and dhp */
    if (fabs(tof - satfrq) >= 0.1) offset(tof,TODEV);
    power(v4,TODEV); 

 status(C);
/* NOE mixing time */
    hsdelay(mix);
 status(D);
/* sampling pulse */
    rgpulse(pw,v1,rof1,rof2);
}
