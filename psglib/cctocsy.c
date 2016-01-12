/* cctocsy - X-nucleus TOCSY with 1H detection				*/
/*  see Fesik et.al., JACS 112,886(1990)                                */
/*									*/
/* pw		- Used for all 90 degree 1H pulses			*/
/* satdly	- Presaturation pulses occur every 10ms for satdly	*/
/* satpwr       - Amplitude for presaturation pulses                    */
/* j            _ X-H coupling constant                                 */
/* pwx		- Used for X-nucleus 90 degree pulses			*/
/* slpwx	- X-nucleus 90-degree pulse width during MLEV-17 period	*/
/* mix		- Total mixing time					*/
/* trim		- trim pulse length					*/
/* tau1		- First two "INEPT-like" delays(1/4j if j<>0)           */
/* tau2		- Refocussing "INEPT-like" delays(1/4j if j<>0)         */
/* pwxpwr	- X-nucleus power during pulses				*/
/* slpwr	- X-nucleus power during MLEVL-17 spin-lock		*/
/* dpwr		- X-nucleus power during acquisition			*/
/* dm		- nny for X-nucleus decoupling during acquisition	*/
/* phase	- 1,2 gives Hypercomplex 2D				*/
/*									*/
/* Phases for sequence are found in tablib/cctocsy			*/
#include <standard.h>

mlev(width,phase1,phase2)
   double width; codeint phase1,phase2;
{
  decrgpulse(width,phase1,0.3e-6,0.0);
  decrgpulse(2*width,phase2,0.3e-6,0.0);
  decrgpulse(width,phase1,0.3e-6,0.0);
}

mlev17(w,x,y,mx,my)
   double w; codeint x,y,mx,my;
{
   mlev(w,my,x); mlev(w,y,mx); mlev(w,y,mx); mlev(w,my,x);
   mlev(w,y,mx); mlev(w,y,mx); mlev(w,my,x); mlev(w,my,x);
   mlev(w,y,mx); mlev(w,my,x); mlev(w,my,x); mlev(w,y,mx);
   mlev(w,my,x); mlev(w,my,x); mlev(w,y,mx); mlev(w,y,mx);
   decrgpulse(2*w,x,0.3e-6,0.0);
}

pulsesequence()
{
   /* Read new parameters */
   double j,pwx,slpwx,slpwr,trim,pwxpwr,mix,tau1,tau2,satdly; int phase;
   pwx=getval("pwx"); slpwx=getval("slpwx"); trim=getval("trim");
   tau1=getval("tau1"); tau2=getval("tau2");
   phase=(int)(getval("phase")+0.5); satdly=getval("satdly"); j=getval("j");
   pwxpwr=getval("pwxpwr"); slpwr=getval("slpwr"); mix=getval("mix");
   if (j!=0.0)
    {
     tau1=1/(4*j); tau2=1/(4*j);
    }

   loadtable("cctocsy");		/* read phase tables		*/

   status(A); 
     initval(pwxpwr,v1); power(v1,DODEV);
     initval(dpwr,v10);
     hsdelay(d1);
   status(B);
     if (satdly!=0.0)			 /* 1H presaturation loop	*/
       { initval(satdly/100,v2);
         starthardloop(v2);
           rgpulse(pw,t1,rof1,0.0);
           delay(.01);
         endhardloop();
       }
     getelem(t2,ct,v3); if (phase==2) incr(v3);
     decrgpulse(pwx,v3,rof1,0.0);
     initval(slpwr,v4); power(v4,DODEV);
     if (d2>0.0)
     {
     delay(d2/2.0);			/* Evolution time		*/
     if (d2 >2.0*pw) rgpulse(2.0*pw,t3,rof1,rof1);
     delay(d2/2.0);
     }
     getelem(t4,ct,v6); add(v6,one,v7); add(v7,one,v8); add(v8,one,v9);
     if (mix>2*66*slpwx)     /* assure twice through hardware loop */
     {
      initval(mix/(66*slpwx),v5);
      decrgpulse(trim,t4,rof1,0.0);
      starthardloop(v5);
       mlev17(slpwx,v6,v7,v8,v9);
      endhardloop();
      decrgpulse(trim,t4,rof1,0.0);
     }
     else
     {
      if (mix>66*slpwx)
       {
        mlev17(slpwx,v6,v7,v8,v9);
       }
      else
       {
        decrgpulse(trim,t4,rof1,0.0);
        decrgpulse(trim,t4,rof1,0.0);
       }
     }
     power(v1,DODEV);
     delay(tau1);
     simpulse(2*pw,2*pwx,t5,t6,rof1,0.0);
     delay(tau1);
     simpulse(pw,pwx,t7,t8,rof1,0.0);
     if (tau2>0.0)
     {
      delay(tau2);
      simpulse(2*pw,2*pwx,t9,t10,rof1,rof2);
      power(v10,DODEV);		/* use dpwr when decoupling*/
      delay(tau2);
     }
   status(C);
     setreceiver(t11);
}
