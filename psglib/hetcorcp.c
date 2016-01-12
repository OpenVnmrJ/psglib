/* hetcorcp - 2D Hetcor for solids

   pw   - 90 degree pulse width
   p1   - tilt pulse (60 degrees recommended)
   trig - ='y'  wait for external trigger
          ='n'  run untriggered
   wim  - # WIM-24 periods used for mixing
*/

#include <standard.h>

pulsesequence()
{
int wim;
char   trig[MAXSTR];
double level1,level1f,level2,level2f;

getstr("trig", trig);
wim=getval("wim");
level1=getval("level1");
level1f=getval("level1f");
level2=getval("level2");
level2f=getval("level2f");

/* set decoupler power levels */
initval(level1,v1);
initval(level1f,v2);  
initval(level2,v3);  
initval(level2f,v4);

initval((double) ix,v9); /* t1 inc */
initval((double) wim,v8);/* # WIM periods */
loadtable("hetcorcp");
setreceiver(t1);

/* start the transient */
status(A);
xmtrphase(zero);
power(v1,DODEV);/* dcplr levels for H-H*/
pwrf(v2,DODEV);
delay(d1);

/*  trigger if required */
if (trig[0] == 'y') xgate(1.0);
rcvroff();

decrgpulse(pw,t2,rof1,0.0);/* prep pulse */

/* loop BB12-BLEW12 ix times  */
/* obs dec obs  dec */
if(ix>0) {
if(ix>1) starthardloop(v9);
simpulse(pw,pw,two,zero,0.0,0.0); /*  x */
simpulse(pw,pw,one,one,0.0,0.0);  /*  y */
simpulse(pw,pw,two,two,0.0,0.0);  /* -x */
simpulse(pw,pw,zero,two,0.0,0.0); /*  y */

simpulse(pw,pw,one,three,0.0,0.0);/*  x */
simpulse(pw,pw,two,zero,0.0,0.0); /*  y */
simpulse(pw,pw,two,two,0.0,0.0);  /* -y */
simpulse(pw,pw,one,one,0.0,0.0);  /* -x */

simpulse(pw,pw,zero,zero,0.0,0.0);/* -y */
simpulse(pw,pw,two,zero,0.0,0.0); /*  x */
simpulse(pw,pw,one,three,0.0,0.0);/* -y */
simpulse(pw,pw,two,two,0.0,0.0);  /* -x */
if(ix>1)endhardloop();
}

decrgpulse(pw,two,0.0,0.0);
decrgpulse(p1,t7,0.0,0.0);

/* WIM 24 Mixing sequence */
if (wim>1) starthardloop(v8);
simpulse(pw,pw,three,t6,0.0,0.0);
simpulse(pw,pw,zero,t3,0.0,0.0);
simpulse(pw,pw,three,t6,0.0,0.0);
simpulse(pw,pw,three,t6,0.0,0.0);
 
simpulse(pw,pw,zero,t3,0.0,0.0);
simpulse(pw,pw,three,t6,0.0,0.0);
simpulse(pw,pw,one,t4,0.0,0.0);
simpulse(pw,pw,zero,t3,0.0,0.0); 
 
simpulse(pw,pw,one,t4,0.0,0.0); 
simpulse(pw,pw,one,t4,0.0,0.0);  
simpulse(pw,pw,zero,t3,0.0,0.0);
simpulse(pw,pw,one,t4,0.0,0.0);
 
simpulse(pw,pw,three,t6,0.0,0.0);
simpulse(pw,pw,two,t5,0.0,0.0);
simpulse(pw,pw,three,t6,0.0,0.0);
simpulse(pw,pw,three,t6,0.0,0.0);

simpulse(pw,pw,two,t5,0.0,0.0);
simpulse(pw,pw,three,t6,0.0,0.0);
simpulse(pw,pw,one,t4,0.0,0.0);
simpulse(pw,pw,two,t5,0.0,0.0);

simpulse(pw,pw,one,t4,0.0,0.0);
simpulse(pw,pw,one,t4,0.0,0.0);  
simpulse(pw,pw,two,t5,0.0,0.0); 
simpulse(pw,pw,one,t5,0.0,0.0);
 
if (wim>1)endhardloop();

status(C);
/* alter decoupler power if required */
if (level1 != level2) power(v3,DODEV);
if (level1f != level2f) pwrf(v4,DODEV);
delay(rof2);
rcvron();
}
