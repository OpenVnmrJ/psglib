/* tocsy1d   -  1D tocsy for rapid scalar coupling correlations;
                internal subtraction of on-resonance and off-
                resonance selective excitation.
                Uses shaped pulse inversion of desired resonance on 
                alternate scans, followed by normal tocsy sequence.

           SELECTIVE PULSE IS DONE VIA PROGR.ATTENUATORS

    ref: a. bax and d.g. davis, j. magn. reson. 65, 355 (1985)

    pw - 90 degree pulse during mlev periods (at power level tpwr)
    p1 - 90 degree excitation pulse (at power p1lvl)
    satdly - presaturation delay using obs. xmtr.
    satfrq - presaturation frequency
    satpwr - presaturation power level
    intsub - 'y': internal subtraction of data acquired by on-resonance
                  and off-resonance selective excitation
             'n': data acquired by on-resonance and off-resonance selec-
                  tive excitation are stored separately
    selfrqon -  frequency of transmitter selective pulse (on resonance);
                array if intsub = 'n'
    selfrqof - frequency of transmitter selective pulse (off resonance);
                an inactive parameter if intsub = 'n'
    selpwr - maximum power of transmitter selective pulse
    shape  -  shape of selective inversion pulse
    seltime - length(in us) of transmitter selective (180) pulse at
		frequency satfrq. (create seltime parameter as a "pulse")
    trim - trim pulse time
    mix - mixing time
    nt - multiple of 32  (intsub = 'n')
         multiple of 64  (intsub = 'y')

            G.Gray  Palo Alto   Sept 1991   */

#include <standard.h>
#define CURRENT 1
#include <shape_pulse.c>

mleva()
{
  txphase(v2); delay(pw);
  txphase(v3); delay(2.0*pw);
  txphase(v2); delay(pw);
}

mlevb()
{
  txphase(v4); delay(pw);
  txphase(v5); delay(2.0*pw);
  txphase(v4); delay(pw);
}


pulsesequence()
{
  /*DEFINE LOCAL VARIABLES */
  double satfrq = getval("satfrq"),
	 satpwr = getval("satpwr"),
	 satdly = getval("satdly"),
	 p1lvl = getval("p1lvl"),
	 trim = getval("trim"),
	 mix = getval("mix"),
	 cycles,
	 selfrqon = getval("selfrqon"),
	 selfrqof,
	 selpwr = getval("selpwr"),
	 seltime = getval("seltime");
  int nisp = 79;
  char intsub[MAXSTR],shape[MAXSTR];
  getstr("shape",shape);
  getstr("intsub",intsub);
  if (intsub[A] == 'y')
    selfrqof = getval("selfrqof");

  /* CALCULATE PHASES */
  if (intsub[A] == 'y')
    hlv(ct,v1);
  else
    assign(ct,v1);
  hlv(v1,v13);
  hlv(v13,v13);
  hlv(v13,v2);
  mod2(v13,v13);
  mod4(v1,v1);
  add(one,v1,v1);
  add(v1,v13,v1);
  assign(v1,oph);
  mod4(v2,v2);
  add(one,v2,v2); 
  add(v13,v2,v2);
  add(one,v2,v3);
  add(one,v3,v4);
  add(one,v4,v5);
  if (intsub[A] == 'y')
  {
    mod2(ct,v14);    /* trigger for the alteration of the saturation freq */
    ifzero(v14);
      add(oph,two,oph);
    endif(v14);
  }
  /* calculate and initialize loop counter */
  cycles = (mix)/(66.0*pw);
  cycles = 2.0*(double)(int)(cycles/2.0);
  mix = cycles*(66.0*pw);
  initval(cycles,v9);

  /* BEGIN ACTUAL PULSE SEQUENCE CODE */
  status(A);
    hsdelay(d1);
  status(B);
    if (satdly>0)
    {
      offset(satfrq,TODEV);
      rlpower(satpwr,TODEV);
      rgpulse(satdly,zero,rof1,rof2); 
      rlpower(p1lvl,TODEV); 
    }

    /* selective pulse or transmitter saturation */
    if (seltime>0)
    {
      if (intsub[A] == 'n')
      {
        offset(selfrqon,TODEV);
      }
      else
      {
        ifzero(v14);
          offset(selfrqon,TODEV);
        elsenz(v14);
          offset(selfrqof,TODEV);
        endif(v14);
      }
      shape_pulse(shape,seltime,zero,t11,t12,selpwr,nisp,rof1,rof1);
      offset(tof,TODEV);
    }
    rlpower(p1lvl,TODEV); 
    rgpulse(p1,v1,30.0e-6,0.2e-6);
    rlpower(tpwr,TODEV); rcvroff(); delay(4e-6);

    /* tocsy mixing time */
    if (mix > 0.0)
    {
      txphase(v5); 
      xmtron(); delay(trim);
      {
        starthardloop(v9);
          mleva(); mlevb(); mlevb(); mleva();
          mlevb(); mlevb(); mleva(); mleva();
          mlevb(); mleva(); mleva(); mlevb();
          mleva(); mleva(); mlevb(); mlevb();
          txphase(v3);
          delay(2*pw);
        endhardloop();
      }
      xmtroff();
    }
    delay(rof2);
    rcvron();
  status(C);
}
