/* SSnoesy -  2D cross relaxation experiment using SS "read" pulse
 
     Noesy for observation of exchangeable protons using a
     Symmetrically-Shifted read pulse.
     Ref. Stephen Smallcombe, JACS 115 4776 (June 2) 1993. 
     Use linear prediction and lsfid to remove lp in multiples of 360 degrees.
     Use "calfa" to set alfa for lp=0 to obtain flat baseline.
     With steering pulses for improved water suppression - 
     Array pwa and pwd for best suppression.
     Positive and negative vaules of pwa and pwd can be used.
     Homospoil in mix can be either normal homospoil, hs='ynn' or
     gradient homospoil, hs='nnn', and gt1>0.
     hsdly hould be set to place homospoil in middle or near end of mix time.
     Axial peaks are displaced via FAD for single transient spectra

*/
#include <standard.h>
pulsesequence()
{
   double mix,SSpwr,pwSS,hsdly,pwa,pwd,gt1,gzlvl1;
   int iphase;
   char sspul[MAXSTR],ptable[MAXSTR],SSshape[MAXSTR];
   getstr("sspul",sspul); getstr("ptable",ptable); getstr("SSshape",SSshape);
   mix=getval("mix"); hsdly=getval("hsdly"); 
   loadtable(ptable);
   iphase = (int) (getval("phase") + 0.5);
   pwa=getval("pwa"); pwd=getval("pwd");
   gt1=getval("gt1"); gzlvl1=getval("gzlvl1");
   SSpwr=getval("SSpwr"); pwSS=getval("pwSS");
     initval(7.0,v2);   /* set up 45 degree phase shift for first pulse    */
     stepsize(45.0,TODEV);
   getelem(t1,ct,v1);     /* 1st pulse phase    */
   getelem(t3,ct,v3);     /* read pulse phase   */
   getelem(t4,ct,v4);     /* receiver phase     */
   if (iphase > 1.5) incr(v1);     /* hypercomplex phase shift */
   assign(v4,oph);
   /*HYPERCOMPLEX MODE USES REDFIELD TRICK TO MOVE AXIAL PEAKS TO EDGE */
    initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v6);
  if ((iphase==1)||(iphase==2)) {add(v1,v6,v1); add(v4,v6,oph);}

    assign(v3,v7);       /* for steering pulses */
    add(v3,one,v8); 
    if (pwa<0.0) { pwa = -pwa; add(v7,two,v7); }
    if (pwd<0.0) { pwd = -pwd; add(v8,two,v8); }

   status(A);
     xmtrphase(v2);
     txphase(v1);
     rlpower(tpwr,TODEV);
     hsdelay(d1);
   status(B);
      rcvroff();
      rgpulse(pw, v1, rof1, 1.0e-6);
      xmtrphase(zero);
      txphase(t2);
      if (d2>0.0)
      delay(d2 -rof1 -(4.0*pw)/3.1414);
      rgpulse(pw, t2, rof1, 1.0e-6);
   status(C);
      txphase(t3);
      rlpower(SSpwr,TODEV);
       delay(hsdly);
    if (gt1>0)
      {
      rgradient('z',gzlvl1);
      delay(gt1);
      rgradient('z',0.0);
      }
       hsdelay(mix-hsdly-gt1);
   status(D);
      shaped_pulse(SSshape,pwSS,t3,rof1,rof1);
      rlpower(SSpwr-12.0,TODEV);
      rgpulse(pwa,v7,2.0e-6,2.0e-6);    /* steering pulses */
      rgpulse(pwd,v8,2.0e-6,2.0e-6);
      rcvron();
}

