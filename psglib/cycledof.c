/* cycledof. Difference NOE experiment

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
              about the frequency "satfrq" at positions given by "spacing"and"pattern" 
              cycle='n' does off-resonance saturation at control
    spacing - spacing(in Hz) of multiplet
    pattern - pattern type ( 1 for singlet, 2 for doublet, etc.)    
    tau     - period spent on a single irradiation point during cycling 
    satpwr - power of decoupler selective pulse
    sattime - length of decoupler selective pulse at frequency satfrq.
    mix - mixing time
    hs - 'ynyn'
    sspul   - sspul='y' does hs-90-hs before d1
    nt - intsub = 'n':  multiple of 16
         intsub = 'y':  multiple of 32


    NOTE:  This pulse sequence requires that the decoupler channel be
           equipped with direct synthesis RF and a linear amplifier.


    s. farmer   5  February    1988
    g. gray  20 june 1988 */


#include <standard.h>

pulsesequence()
{
/*DEFINE LOCAL VARIABLES */
  double mix,satfrq,control,satpwr,sattime,
         tau,spacing;
  int pattern,times,jj;
  char intsub[MAXSTR],cycle[MAXSTR],sspul[MAXSTR];


/* LOAD AND INITIALIZE VARIABLES */
  tau = getval("tau");
  getstr("intsub",intsub); getstr("cycle",cycle); getstr("sspul",sspul); 
  satfrq = getval("satfrq"); control = getval("control");
  sattime = getval("sattime"); spacing = getval("spacing");
  satpwr = getval("satpwr"); pattern = (int)getval("pattern");
  mix = getval("mix");
  if (pattern == 0) pattern = 1; if (tau == 0.0) tau = 0.1;
  times = (int)(sattime/(pattern*tau));
  initval(dhp,v7);
  initval(satpwr,v6);
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
    hsdelay(hst); pulse(pw,zero); hsdelay(hst);
    }
    hsdelay(d1);
    power(v5,TODEV);  delay(0.2e-6); /*reduce leakage via xmtr*/
 status(B);

/* selective pulse or decoupler saturation */
    power(v6,DODEV);
    /* no cycling or interleaved subtraction (make control an array)*/
    if ((intsub[0] == 'n') && (cycle[0] == 'n'))
     {
      offset(satfrq,DODEV);
      delay(sattime);
     }
    /* interleaved subtraction without cycling */
    if ((intsub[0] == 'y') && (cycle[0] == 'n'))
      {
       ifzero(v14);
          offset(control,DODEV);
          delay(sattime);
       elsenz(v14);
          offset(satfrq,DODEV);
          delay(sattime);
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
         offset(startfrq,DODEV);
         delay(tau);
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
           offset(startfrq,DODEV);
           delay(tau);
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
           offset(startfrq,DODEV);
           delay(tau);
           startfrq = startfrq + spacing;
          }
        }  
      endif(v14);
     }
    /* restore power levels as controlled by tpwr and dhp */
    offset(dof,DODEV);
    power(v7,DODEV); 
    power(v4,TODEV); delay(5.0e-5);

 status(C);
/* NOE mixing time */
    hsdelay(mix);
 status(D);
/* sampling pulse */
    rgpulse(pw,v1,rof1,rof2);
}
