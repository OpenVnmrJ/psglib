/* gmqcosyps
	phase-sensitive pulsed gradient enhanced mq cosy
	(currently (930326) this is set up with the gmqcosy macro)
   
   Parameters:
	gzlvl1 - try 15000
	gt1 - try 0.0025 seconds
	phase - 1,2 (run and process as a standard hypercomplex dataset)

   Processing:
	if phase=1,2: use wft2da
		(after running the gmqcosy macro to setup, type rp1=90)
	 
	pak 920506 -
	eh 9301 - tested and working (was ehgmqcosy)
	pak 930326 - added more documentation
*/


#include <standard.h>

pulsesequence()
{
        double gzlvl1,qlvl,grise,gstab,gt1,satpwr,satdly,phase;
	int icosel,
            iphase;
        satpwr = getval("satpwr");
   	satdly = getval("satdly");
        grise=getval("grise");
        gstab=getval("gstab");
        gt1=getval("gt1");
        gzlvl1=getval("gzlvl1");
        qlvl=getval("qlvl");
        phase=getval("phase");
        iphase = (int) (phase + 0.5);


/* DETERMINE STEADY-STATE MODE */
   assign(oph,v2);
   assign(oph,v1);
   if (iphase == 2)
      incr(v1);

/* HYPERCOMPLEX MODE */
   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v6);
   if ((iphase==1) || (iphase==2))
      {add(v1,v6,v1); add(oph,v6,oph);} 
	

     status(A);
	rlpower(satpwr,TODEV);
	delay(d1);
	rgpulse(satdly,zero,rof1,rof2);
	rlpower(tpwr,TODEV);
        rcvroff();
	rgpulse(pw,v1,rof1,rof2);
     status(B);
        if (d2 > rof1 + 4.0*pw/3.1416)
           delay(d2 - rof1 - 4.0*pw/3.1416);
     status(C);
        rgpulse(pw,v2,rof1,rof2);
      
        delay(gt1+2.0*grise+24.4e-6);
        rgpulse(2.0*pw,v2,rof1,rof2);
        icosel=-1;
        rgradient('z',gzlvl1*(double)icosel);
	delay(gt1+grise);
	rgradient('z',0.0);
        txphase(oph);
	delay(grise);

	rgpulse(pw,v2,rof1,rof2);

        rgradient('z',gzlvl1*qlvl);
	delay(gt1+grise);
	rgradient('z',0.0);
	delay(grise);
	delay(gstab);
     status(D);
}
