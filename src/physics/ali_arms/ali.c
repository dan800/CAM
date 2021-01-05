#include <stdio.h>
#include <stdlib.h>
#include <./include/struct.h>
#include <./include/constants.h>
#include <./include/subroutines.h>
#include <./include/functions.h>
#include <./include/mv_utils.h>


int ali_(double *H,
         double *P,
         double *T,
         double *CO2,
         double *O,
         double *N2,
         double *O2,
         double *CH_SUM,
         int    *N)

{
 PARAMETERINFO 	*pars;
 PARAMETERINFO	parsdata;
 DEBYE_INFO	*deb;

 MOLECULE	*mol;
 VIBLEVEL	*v_level;

 BAND_CVT	*cvt_band;
 PBAND		pcvt_band;

 BAND_CVV	*cvv_band;

 ATMOSPHERE	atmos;
 INTEGR_ODF 	integr_ODF;

 BAND_ODF	*band_odf;
 
 LINE_ODF	*line_odf;

 double		**PopV, **PopVN, **DeltaV, **Qrot;
 double		**Rmat, **RmatBig;
 double		***PopStack, ***PopStackDelta;
 
 double		time;

 double		t0,t1,t2,t3,t4,t5;
 
 double		maxDelta,maxDeltaPrev;

 int		NMOL,ND,NVL,NCVT,NCVV,N15MKM;
 int		id,ik,iter;
 int		KBACK,KBACK1;

 PopV=NULL;PopVN=NULL;DeltaV=NULL;Qrot=NULL;
 Rmat=NULL;RmatBig=NULL;PopStack=NULL;PopStackDelta=NULL;

 pars=NULL;deb=NULL;mol=NULL;v_level=NULL;
 cvt_band=NULL;cvv_band=NULL;
 
 int ii;

 pars=&parsdata;
 memset(pars,0,sizeof(parsdata));
 
 init_parameters(pars);
     
#ifdef PRINT 
 printf("\nThe program starts with the following parameters:\n");
 print_PARAMETERINFO(pars);
#endif
 
/*    
 * mol=init_molecule_data(mol,pars);
 * NMOL=pars->NMOL;
 */
    
    ND=N[0];
    
    for(ii=0;ii<ND-1;ii++){
        CH_SUM[ii] = -0.0012;
    }
    
    ii = 20;
    CH_SUM[ii] = H[ii];
    ii = 10;
    CH_SUM[ii] = CO2[ii];
    
    return 0;
}
