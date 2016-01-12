/*
    superb-w pulse sequences  [Fetler et al., JMR Series B, 101, 17, 1993;
    Fetler et al., Biophys. J., 68, 681, 1995]

    The SUPERB-W pulse sequences are composite solvent suppression 
    pulse sequences, and are similar to the binomial pulse sequences
    except that they are self-refocusing in the excited bandwidth.  
    They are easy to use.  The following parameters need to be specified,
    and may be created with the superw macro:

	iseq - selects the pulse sequence to be used
	fmax - sets the bandwidth of uniform phase excitation
	tpwr - power level 
	pws  - pulse length (use at least 25 usec)
	rof1 - pulse switching time (use a short value, i.e. 2 usec,
		to obtain a smaller lp)

    Most pulses are given as a fraction times pws (0.xxxx*pws), and the
    overall flip angle is given by pws.  tpwr should be set appropriately 
    to obtain a given flip angle.  Short pulses less than 0.01*pws
    are commented out, and pws > 25 usec is recommended. 

    The value of iseq generally consists of the first digit, M,
    specifying the null order plus one, and the last two digits, N,
    specifying the "number" of pulses in the pulse sequence.  

	iseq = 100*M + N                               (1)

    Allowed values of iseq:  121, 205, 308, 311, 416, 423, 432
    [121 from Starcuk and Sklenar, JMR, 61, 567, 1985;
    this is an exception to eq. (1)].  

    The region of uniform phase excitation covers from approximately 
    M*f0/N to (2*N-M)*f0/N, where f0 is the center of uniform excitation 
    region, in hertz from solvent.  As with the binomial pulse sequences, 
    the delay "tau" is equal to 1/(2*f0).  The delay tau is set by fmax
    as follows:

	fmax < 0: f0 set to sw/4
	fmax > 0: fmax set to (2*N-M)*f0/N

    Automatic calculation of the null and excitation bandwidth may be
    done using the macro superw('calc') or superwcalc in maclib.  

    For typical proton work, try fmax = -100 or fmax = 3000, with tof
    on resonance with solvent, and iseq = 311 or iseq = 432.  You may
    wish to use solvent subtraction data processing to further reduce
    the solvent signal (ssfilter, ssntaps).

*/

#include <standard.h>
pulsesequence()
{
   double tau, fmax, pws;
   int iseq, null, pslen;
   iseq = getval("iseq"); fmax = getval("fmax");
   pws = getval("pws"); 

/* Calculate tau */
   if (fmax <= 0.1)
     {
       tau = 2.0/getval("sw");
     }
   else
     {
     if (iseq == 121)
       {
         tau = 1.0/(2.0*fmax);
       }
     else
       {
       /* iseq = 100*null + pslen */
       null = iseq / 100;
       pslen = iseq % 100;
       tau = (double) (2*pslen - null)/((double) (2*pslen) * fmax);
       }
     }

   assign(oph,v8); add(v8,two,v9);

   status(A);
     hsdelay(d1);

   status(B);

  /*
  Pulse sequence is (1,-2,1)
  Number of "pulses" is    3
  */

  if (iseq == 121)
     {
     rgpulse(  0.2500*pws,v8,rof1,rof1);
       delay(tau- 3*rof1-  0.3750*pws);
     rgpulse(  0.2500*pws,v9,rof1,rof1);
     rgpulse(  0.2500*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.3750*pws);
     rgpulse(  0.2500*pws,v8,rof1,rof1);
     }


  /*
  Pulse sequence is SUPERB-W(0, 2pi/5)
  Number of "pulses" is    5
  */

  if (iseq == 205)
     {
     rgpulse(  0.1910*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.2500*pws);
     rgpulse(  0.3090*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.3090*pws);
     rgpulse(  0.3090*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.2500*pws);
     rgpulse(  0.1910*pws,v9,rof1,rof1);
       delay(tau- 5*rof1-  0.5955*pws);
     rgpulse(  0.1910*pws,v8,rof1,rof1);
     rgpulse(  0.3090*pws,v8,rof1,rof1);
     rgpulse(  0.3090*pws,v8,rof1,rof1);
     rgpulse(  0.1910*pws,v8,rof1,rof1);
     }


  /*
  Pulse sequence is SUPERB-W(1, 3pi/8)
  Number of "pulses" is    8
  */
 
  if (iseq == 308)
     {
     rgpulse(  0.0864*pws,v8,rof1,rof1);
     rgpulse(  0.1221*pws,v8,rof1,rof1);
       delay(tau- 4*rof1-  0.2358*pws);
     rgpulse(  0.2328*pws,v8,rof1,rof1);
     rgpulse(  0.0303*pws,v8,rof1,rof1);
       delay(tau- 3*rof1-  0.1841*pws);
     rgpulse(  0.1050*pws,v8,rof1,rof1);
       delay(tau- 3*rof1-  0.1482*pws);
     rgpulse(  0.1050*pws,v9,rof1,rof1);
     rgpulse(  0.0864*pws,v9,rof1,rof1);
       delay(tau- 4*rof1-  0.3249*pws);
     rgpulse(  0.0303*pws,v9,rof1,rof1);
     rgpulse(  0.4282*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.4961*pws);
     rgpulse(  0.5337*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.4444*pws);
     rgpulse(  0.1221*pws,v9,rof1,rof1);
     rgpulse(  0.2328*pws,v9,rof1,rof1);
       delay(tau- 4*rof1-  0.6585*pws);
     rgpulse(  0.4282*pws,v8,rof1,rof1);
     rgpulse(  0.5337*pws,v8,rof1,rof1);
     }

 
  /*
  Pulse sequence is SUPERB-W(1, 3pi/11)
  Number of "pulses" is   11
  */
 
  if (iseq == 311)
     {
     rgpulse(  0.1089*pws,v8,rof1,rof1);
       delay(tau- 3*rof1-  0.1416*pws);
     rgpulse(  0.0852*pws,v8,rof1,rof1);
     rgpulse(  0.0891*pws,v8,rof1,rof1);
       delay(tau- 4*rof1-  0.1709*pws);
     rgpulse(  0.0456*pws,v8,rof1,rof1);
     rgpulse(  0.1220*pws,v8,rof1,rof1);
       delay(tau- 4*rof1-  0.1257*pws);
     rgpulse(  0.0549*pws,v8,rof1,rof1);
     rgpulse(  0.0290*pws,v8,rof1,rof1);
       delay(tau- 3*rof1-  0.0694*pws);
     rgpulse(  0.0549*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.1310*pws);
     rgpulse(  0.1220*pws,v9,rof1,rof1);
     rgpulse(  0.0852*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.2659*pws);
     rgpulse(  0.3247*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.3462*pws);
     rgpulse(  0.3676*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.3425*pws);
     rgpulse(  0.0891*pws,v9,rof1,rof1);
     rgpulse(  0.2283*pws,v9,rof1,rof1);
       delay(tau- 5*rof1-  0.2504*pws);
     rgpulse(  0.0290*pws,v9,rof1,rof1);
     rgpulse(  0.1089*pws,v9,rof1,rof1);
     rgpulse(  0.0456*pws,v9,rof1,rof1);
       delay(tau- 6*rof1-  0.5520*pws);
     rgpulse(  0.2283*pws,v8,rof1,rof1);
     rgpulse(  0.3247*pws,v8,rof1,rof1);
     rgpulse(  0.3676*pws,v8,rof1,rof1);
     }

 
  /*
  Pulse sequence is SUPERB-W(2, 4pi/16)
  Number of "pulses" is   16
  */
 
  if (iseq == 416)
     {
     rgpulse(  0.0865*pws,v9,rof1,rof1);
     rgpulse(  0.0147*pws,v9,rof1,rof1);
       delay(tau- 4*rof1-  0.1281*pws);
     rgpulse(  0.0542*pws,v9,rof1,rof1);
     rgpulse(  0.1007*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.1432*pws);
     rgpulse(  0.1315*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.0813*pws);
/*     rgpulse(  0.0019*pws,v9,rof1,rof1); */
     rgpulse(  0.0292*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.0733*pws);
     rgpulse(  0.0292*pws,v8,rof1,rof1);
     rgpulse(  0.0865*pws,v8,rof1,rof1);
       delay(tau- 4*rof1-  0.1863*pws);
     rgpulse(  0.1007*pws,v8,rof1,rof1);
     rgpulse(  0.1563*pws,v8,rof1,rof1);
       delay(tau- 3*rof1-  0.2978*pws);
     rgpulse(  0.3385*pws,v8,rof1,rof1);
       delay(tau- 4*rof1-  0.3303*pws);
     rgpulse(  0.0614*pws,v8,rof1,rof1);
     rgpulse(  0.2393*pws,v8,rof1,rof1);
     rgpulse(  0.0213*pws,v8,rof1,rof1);
       delay(tau- 6*rof1-  0.2612*pws);
     rgpulse(  0.0147*pws,v8,rof1,rof1);
     rgpulse(  0.1315*pws,v8,rof1,rof1);
     rgpulse(  0.0542*pws,v8,rof1,rof1);
       delay(tau- 3*rof1-  0.1012*pws);
/*     rgpulse(  0.0019*pws,v8,rof1,rof1); */
       delay(tau- 2*rof1-  0.1098*pws);
     rgpulse(  0.1563*pws,v9,rof1,rof1);
     rgpulse(  0.0614*pws,v9,rof1,rof1);
       delay(tau- 5*rof1-  0.3048*pws);
     rgpulse(  0.0213*pws,v9,rof1,rof1);
     rgpulse(  0.3385*pws,v9,rof1,rof1);
     rgpulse(  0.0321*pws,v9,rof1,rof1);
       delay(tau- 4*rof1-  0.4283*pws);
     rgpulse(  0.4648*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.4374*pws);
     rgpulse(  0.4100*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.3247*pws);
     rgpulse(  0.2393*pws,v9,rof1,rof1);
       delay(tau- 4*rof1-  0.5731*pws);
     rgpulse(  0.0321*pws,v8,rof1,rof1);
     rgpulse(  0.4100*pws,v8,rof1,rof1);
     rgpulse(  0.4648*pws,v8,rof1,rof1);
     }

 
  /*
  Pulse sequence is SUPERB-W(2, 4pi/23)
  Number of "pulses" is   23
  */
 
  if (iseq == 423)
     {
     rgpulse(  0.0381*pws,v9,rof1,rof1);
     rgpulse(  0.0159*pws,v9,rof1,rof1);
       delay(tau- 4*rof1-  0.0748*pws);
     rgpulse(  0.0605*pws,v9,rof1,rof1);
     rgpulse(  0.0351*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.1051*pws);
/*     rgpulse(  0.0089*pws,v9,rof1,rof1); */
     rgpulse(  0.1056*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.1098*pws);
     rgpulse(  0.1051*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.0858*pws);
     rgpulse(  0.0561*pws,v9,rof1,rof1);
     rgpulse(  0.0105*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.0354*pws);
/*     rgpulse(  0.0042*pws,v9,rof1,rof1); */
       delay(tau- 2*rof1-  0.0381*pws);
     rgpulse(  0.0159*pws,v8,rof1,rof1);
     rgpulse(  0.0561*pws,v8,rof1,rof1);
       delay(tau- 4*rof1-  0.1106*pws);
     rgpulse(  0.0351*pws,v8,rof1,rof1);
     rgpulse(  0.1051*pws,v8,rof1,rof1);
/*     rgpulse(  0.0089*pws,v8,rof1,rof1); */
       delay(tau- 4*rof1-  0.1808*pws);
     rgpulse(  0.1056*pws,v8,rof1,rof1);
     rgpulse(  0.1069*pws,v8,rof1,rof1);
       delay(tau- 4*rof1-  0.2312*pws);
     rgpulse(  0.2226*pws,v8,rof1,rof1);
     rgpulse(  0.0272*pws,v8,rof1,rof1);
       delay(tau- 3*rof1-  0.2508*pws);
     rgpulse(  0.2462*pws,v8,rof1,rof1);
/*     rgpulse(  0.0056*pws,v8,rof1,rof1); */
       delay(tau- 4*rof1-  0.2335*pws);
     rgpulse(  0.0199*pws,v8,rof1,rof1);
     rgpulse(  0.1822*pws,v8,rof1,rof1);
     rgpulse(  0.0131*pws,v8,rof1,rof1);
       delay(tau- 6*rof1-  0.1787*pws);
     rgpulse(  0.0105*pws,v8,rof1,rof1);
     rgpulse(  0.0714*pws,v8,rof1,rof1);
     rgpulse(  0.0605*pws,v8,rof1,rof1);
       delay(tau- 4*rof1-  0.0923*pws);
/*     rgpulse(  0.0042*pws,v8,rof1,rof1); */
     rgpulse(  0.0381*pws,v8,rof1,rof1);
       delay(tau- 2*rof1-  0.0568*pws);
     rgpulse(  0.0714*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.1268*pws);
     rgpulse(  0.1822*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.2278*pws);
     rgpulse(  0.0272*pws,v9,rof1,rof1);
     rgpulse(  0.2462*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.3021*pws);
     rgpulse(  0.3307*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.3377*pws);
     rgpulse(  0.3446*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.3282*pws);
/*     rgpulse(  0.0056*pws,v9,rof1,rof1); */
     rgpulse(  0.3061*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.2737*pws);
     rgpulse(  0.0131*pws,v9,rof1,rof1);
     rgpulse(  0.2226*pws,v9,rof1,rof1);
       delay(tau- 4*rof1-  0.1813*pws);
     rgpulse(  0.1069*pws,v9,rof1,rof1);
     rgpulse(  0.0199*pws,v9,rof1,rof1);
       delay(tau- 5*rof1-  0.5541*pws);
     rgpulse(  0.3061*pws,v8,rof1,rof1);
     rgpulse(  0.3307*pws,v8,rof1,rof1);
     rgpulse(  0.3446*pws,v8,rof1,rof1);
     }

 
  /*
  Pulse sequence is SUPERB-W(2, 4pi/32)
  Number of "pulses" is   32
  */
 
  if (iseq == 432)
     {
/*     rgpulse(  0.0072*pws,v9,rof1,rof1); */
     rgpulse(  0.0207*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.0403*pws);
     rgpulse(  0.0514*pws,v9,rof1,rof1);
/*     rgpulse(  0.0012*pws,v9,rof1,rof1); */
       delay(tau- 3*rof1-  0.0617*pws);
     rgpulse(  0.0129*pws,v9,rof1,rof1);
     rgpulse(  0.0579*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.0756*pws);
     rgpulse(  0.0803*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.0799*pws);
     rgpulse(  0.0780*pws,v9,rof1,rof1);
/*     rgpulse(  0.0015*pws,v9,rof1,rof1); */
       delay(tau- 2*rof1-  0.0736*pws);
     rgpulse(  0.0677*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.0567*pws);
     rgpulse(  0.0350*pws,v9,rof1,rof1);
     rgpulse(  0.0106*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.0303*pws);
     rgpulse(  0.0149*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.0185*pws);
/*     rgpulse(  0.0000*pws,v8,rof1,rof1); */
     rgpulse(  0.0149*pws,v8,rof1,rof1);
/*     rgpulse(  0.0072*pws,v8,rof1,rof1); */
       delay(tau- 3*rof1-  0.0421*pws);
     rgpulse(  0.0106*pws,v8,rof1,rof1);
     rgpulse(  0.0514*pws,v8,rof1,rof1);
       delay(tau- 4*rof1-  0.0816*pws);
/*     rgpulse(  0.0012*pws,v8,rof1,rof1); */
     rgpulse(  0.0586*pws,v8,rof1,rof1);
     rgpulse(  0.0414*pws,v8,rof1,rof1);
       delay(tau- 4*rof1-  0.1186*pws);
     rgpulse(  0.0579*pws,v8,rof1,rof1);
     rgpulse(  0.0780*pws,v8,rof1,rof1);
       delay(tau- 4*rof1-  0.1493*pws);
     rgpulse(  0.0541*pws,v8,rof1,rof1);
     rgpulse(  0.1085*pws,v8,rof1,rof1);
       delay(tau- 4*rof1-  0.1704*pws);
     rgpulse(  0.0151*pws,v8,rof1,rof1);
     rgpulse(  0.1631*pws,v8,rof1,rof1);
       delay(tau- 4*rof1-  0.1795*pws);
     rgpulse(  0.0145*pws,v8,rof1,rof1);
     rgpulse(  0.1662*pws,v8,rof1,rof1);
       delay(tau- 4*rof1-  0.1749*pws);
     rgpulse(  0.0177*pws,v8,rof1,rof1);
     rgpulse(  0.1514*pws,v8,rof1,rof1);
       delay(tau- 4*rof1-  0.1562*pws);
/*     rgpulse(  0.0015*pws,v8,rof1,rof1); */
     rgpulse(  0.0803*pws,v8,rof1,rof1);
     rgpulse(  0.0615*pws,v8,rof1,rof1);
       delay(tau- 5*rof1-  0.1240*pws);
     rgpulse(  0.0241*pws,v8,rof1,rof1);
     rgpulse(  0.0677*pws,v8,rof1,rof1);
     rgpulse(  0.0129*pws,v8,rof1,rof1);
       delay(tau- 5*rof1-  0.0802*pws);
     rgpulse(  0.0207*pws,v8,rof1,rof1);
     rgpulse(  0.0350*pws,v8,rof1,rof1);
       delay(tau- 2*rof1-  0.0279*pws);
/*     rgpulse(  0.0000*pws,v9,rof1,rof1); */
       delay(tau- 1*rof1-  0.0293*pws);
     rgpulse(  0.0586*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.0871*pws);
     rgpulse(  0.0615*pws,v9,rof1,rof1);
     rgpulse(  0.0541*pws,v9,rof1,rof1);
       delay(tau- 4*rof1-  0.1410*pws);
     rgpulse(  0.1514*pws,v9,rof1,rof1);
     rgpulse(  0.0151*pws,v9,rof1,rof1);
       delay(tau- 4*rof1-  0.1868*pws);
     rgpulse(  0.1662*pws,v9,rof1,rof1);
     rgpulse(  0.0409*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.2207*pws);
     rgpulse(  0.2342*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.2398*pws);
     rgpulse(  0.2454*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.2423*pws);
     rgpulse(  0.2393*pws,v9,rof1,rof1);
       delay(tau- 2*rof1-  0.2278*pws);
     rgpulse(  0.2163*pws,v9,rof1,rof1);
       delay(tau- 3*rof1-  0.1970*pws);
     rgpulse(  0.1631*pws,v9,rof1,rof1);
     rgpulse(  0.0145*pws,v9,rof1,rof1);
       delay(tau- 4*rof1-  0.1519*pws);
     rgpulse(  0.1085*pws,v9,rof1,rof1);
     rgpulse(  0.0177*pws,v9,rof1,rof1);
       delay(tau- 4*rof1-  0.0958*pws);
     rgpulse(  0.0414*pws,v9,rof1,rof1);
     rgpulse(  0.0241*pws,v9,rof1,rof1);
       delay(tau- 7*rof1-  0.5208*pws);
     rgpulse(  0.0409*pws,v8,rof1,rof1);
     rgpulse(  0.2163*pws,v8,rof1,rof1);
     rgpulse(  0.2342*pws,v8,rof1,rof1);
     rgpulse(  0.2393*pws,v8,rof1,rof1);
     rgpulse(  0.2454*pws,v8,rof1,rof1);
     }
 
}
