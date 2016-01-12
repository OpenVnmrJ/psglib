/* 	CIGAR - Improved IMPEACH
		Introduces a Varaible J Scale for JHH skewing

        Features included:
                F1 Axial Displacement
                Randomization of Magnetization prior to relaxation delay
                        with G-90-G
                        [selected by sspul flag]
                Includes a 2-stage J filter and BIRD based one-bond correlation
                 	suppression

        Paramters:
                sspul :         y - selects magnetization randomization option
                hsglvl:         Homospoil gradient level (DAC units)
                hsgt    :       Homospoil gradient time
                gzlvl1  :       encoding Gradient level
                gt1     :       encoding gradient time
                gzlvl3  :       decoding Gradient level
                gt3     :       decoding gradient time
                gzlvl0  :	gradinets used for echos
                gt0	:	time for echo gradients
                gstab   :       recovery delay
        j1max & j1min   :       One-bond XH coupling constants
                jnxh    :       multiple bond XH coupling constant for gHMBC
        jnmax & jnmin	:	multiple bond XH coupling constants for accordion
        	jscaleU :	scaling factor for JHH skewing [typically 0]
                pwxlvl  :       X-nucleus pulse power
                pwx     :       X-nucleus 90 deg pulse width
                d1      :       relaxation delay
                d2      :       Evolution delay

KrishK	-	CIGAR-HMBC		: February '99
KrishK	-	can now be used for "static" gHMBC as well.  : Feb. 00
				
*/


#include <standard.h>
static double d2_init = 0.0;

static int ph1[1] = {0};
static int ph3[4] = {0,2};
static int ph4[1] = {0};
static int ph5[8] = {0,0,2,2};
static int ph6[8] = {0,2,2,0};

pulsesequence()
{
  double j1min,
  	 j1max,
         pwxlvl,
         pwx,
         gzlvl0,
         gt0,
	 gzlvl1,
	 gt1,
         gzlvl3,
         gt3,
	 bird,
         gstab,
         hsglvl,
         hsgt,
         tauX,
         tau,
         tautau,
         t1max,
	 tauA,
	 tauB,
         taumb;
  char   sspul[MAXSTR],
  	 accord[MAXSTR];
  int    ijscaleU,
         t1_counter;

  j1min = getval("j1min");
  j1max = getval("j1max");
  pwxlvl = getval("pwxlvl");
  pwx = getval("pwx");
  getstr("sspul",sspul);
  getstr("accord",accord);
  gzlvl0 = getval("gzlvl0");
  gt0 = getval("gt0");
  gzlvl1 = getval("gzlvl1");
  gt1 = getval("gt1");
  gzlvl3 = getval("gzlvl3");
  gt3 = getval("gt3");
  gstab = getval("gstab");
  hsglvl = getval("hsglvl");
  hsgt = getval("hsgt");
  tauX = 1/(2*(getval("jnmax")));
  tauA = 1/(2*(j1min + 0.146*(j1max - j1min)));
  tauB = 1/(2*(j1max - 0.146*(j1max - j1min)));
  bird = (tauA+tauB)/2;
  ijscaleU = (int)(getval("jscaleU") + 0.5);
  t1max = ((getval("ni")-1)/getval("sw1")) + 0.0005;
  tautau = ((ijscaleU - 1)*d2);
  taumb = 1/(2*(getval("jnmin")));
   if(ix == 1)
     d2_init = d2;
   t1_counter = (int) ( (d2-d2_init)*sw1 + 0.5);

  if (accord[0] == 'y')
  {
  if ((taumb - tauX) < t1max/2)
    taumb = tauX + t1max/2;  
  tau = ((taumb - tauX)*t1_counter/ni);
  }
  else
  {
    taumb = 1/(2*(getval("jnxh")));
    if (ijscaleU > 0) 
      tau = 0.0;
    else
      tau = t1max/2;
  }
  
  settable(t1,1,ph1);
  settable(t3,2,ph3);
  settable(t4,1,ph4);
  settable(t5,4,ph5);
  settable(t6,4,ph6);

  getelem(t3,ct,v3);
  getelem(t6,ct,oph);

  initval(2.0*(double)((int)(d2*getval("sw1")+0.5)%2),v10);
  add(v3,v10,v3);
  add(oph,v10,oph);

  status(A);
     decpower(pwxlvl);
     if (sspul[0] == 'y')
     {
        zgradpulse(hsglvl,hsgt);
        rgpulse(pw,zero,rof1,rof1);
        zgradpulse(hsglvl,hsgt);
     }

     delay(d1);
     rcvroff();

  status(B);
     rgpulse(pw,t1,rof1,rof2);

/* Start of J filter  */
     zgradpulse(gzlvl0/2,gt0);
     delay(tauA - gt0);
     decrgpulse(pwx, zero, rof1, rof1);
     zgradpulse(-gzlvl0/3,gt0);
     delay(tauB - gt0);
     decrgpulse(pwx, zero, rof1, rof1);
     zgradpulse(-gzlvl0/6,gt0);
     delay(gstab);
/* End of J filter */

/* Start of BIRD suppression  */

        zgradpulse(gzlvl0*1.5,gt0);
        delay(gstab);
        rgpulse(pw,zero,rof1,rof1);
        delay(bird);
        simpulse(2*pw,2*pwx,zero,zero,rof1,rof1);
        delay(bird);
        rgpulse(pw,zero,rof1,rof2);
        zgradpulse(-gzlvl0*1.5,gt0);
        delay(gstab);
/* End of BIRD suppression */

/*    Start of CT-VD */

    if (!((accord[0] == 'n') && (ijscaleU == 1)))
      {
       delay(tau/2 + tautau/4);
        zgradpulse(gzlvl0,gt0);
        delay(gstab);
        decrgpulse(2*pwx,two,rof1,rof1);
        zgradpulse(-gzlvl0,gt0);
        delay(gstab);
       delay(tau/2 + tautau/4);
      }
       
    if (accord[0] == 'y')
     delay(taumb - tau);
    else
     delay(taumb);
     
/* End of CT-VD */

     decrgpulse(pwx,v3,rof1,rof1);
     zgradpulse(gzlvl1,gt1);
     delay(gstab);
     delay(d2/2);
     rgpulse(pw*2.0,t4,rof1,rof1);
     delay(d2/2);
     zgradpulse(gzlvl1,gt1);
     delay(gstab);
     decrgpulse(pwx,t5,rof1,rof1);

     zgradpulse(gzlvl3,gt3);
     decpower(dpwr);
     delay(t1max/2 + tautau/2);
     rcvron();
  status(C);
}