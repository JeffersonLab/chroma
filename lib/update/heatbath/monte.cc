// $Id: monte.cc,v 3.0 2006-04-03 04:59:07 edwards Exp $
/*! \file
 *  \brief Routine for doing the hybrid (monte carlo) algorithm. 
 */

#error "Not tested (or even compiled). However, reasonably well converted."


#include "chromabase.h"
#include "update/heatbath/mciter32.h"
#include "meas/glue/mesplq.h"
#include "meas/glue/fuzglue.h"
#include "meas/glue/wilslp.h"
#include "meas/glue/topol.h"
#include "meas/smear/smear.h"
#include "meas/gfix/gfix.h"

namespace Chroma {

//! Routine for doing the hybrid (monte carlo) algorithm.
/*!
 * \ingroup heatbath
 *
 *
 * \param u           gauge field ( Modify )
 * \param NOver       number of overrelaxations steps ( Read )
 * \param NHeat       number of heatbath trials ( Read )
 *
 * \param j_decay     direction of exponential decay --- for sideways potl ! ( Read )
 *
 * \param GlueP       flag for glueball measurements ( Read )
 * \param WilsonP     flag for Wilson loop measurements ( Read )
 * \param TopolP      flag for topology measurements ( Read )
 * \param PrtMes      measurements per printout to RESULT file
 * \param MesItr      number of iterations in mod for measurements ( Read )
 * \param NumTry      number of heatbath link trials performed ( Write )
 * \param NumFail     number of heatbath link trials failed ( Write )
 * \param NumItr      number of iterations to perform ( Read ) 
 * \param TotalItr    total number of iterations to perform ( Read ) 
 * \param BlkMax      maximum number of SU(3) trace maximimations allowed ( Read )
 * \param BlkAccu     blocking accuracy in glueball code ( Read )
 * \param TopolItr    number of iterations in mod for topology measurement ( Read )
 * \param TopAccu     accuracy for convergence of topological charge ( Read )
 * \param ActAccu     accuracy for convergence of action ratio ( Read )
 * \param NumTop      number of topological charge measurements ( Read )
 * \param NumCool     number of cooling sweeps per topological charge measurement ( Read )
 * \param HQPotP      flag for heavy quark potential measurements ( Read )
 * \param HQPotItr    number of iterations in mod for potential measurement ( Read )
 * \param FacSmea     smearing factor for potential measurement ( Read )
 * \param NumSmea     number of smearing iterations in potential measurement ( Read )
 * \param GFixP       flag for Landau gauge fixing ( Read )
 * \param GFixItr     number of iterations in mod for Landau gauge fixing ( Read )
 * \param ORlxDo      do over-relaxation in gauge fixing ( Read )
 * \param OrPara      over-relaxation parameter in gauge-fixing ( Read ) 
 */

void monte(multi1d<LatticeColorMatrix>& u(Nd),	/* New gauge field cfg. */
	   int NOver,			/* Number of overrelaxations steps */
	   int NHeat,			/* Number of heatbath trials */
	   int j_decay,		/* Direction to measure sideways potl */
                                /* global t_dir is used for everything else! */
	   int GlueP,			/* Measure glueballs on-line? */
	   int WilsonP,		/* Measure Wilson loops on-line? */
	   int& NumTry,		/* Number of heatbath link trials for this run */
	   int& NumFail,	/* Number of heatbath link failures for this run */
	   int PrtMes,		/* Measurements per printout to RESULT file */
	   int MesItr,		/* Iterations per measurement */
	   int NumItr,		/* Number of iterations for this run */
	   int TotalItr,		/* Total number of iterations */

	   int TopolP,		/* Measure topology on-line? */
	   int TopolItr,		/* Iterations per topology measurement */
	   Real TopAccu,                  /* Topological charge measurement accuracy */
	   Real ActAccu,                  /* Action measurement accuracy */
	   int NumTop,		/* Frequency of charge measurement */
	   int NumCool,		/* Total number of cooling sweeps */

	   Real BlkAccu,		        /* Accuracy of blocking in fuzglue. */
	   int BlkMax,		/* Maximum number of proj. iterations in fuzglue. */

	   int HQPotP,		/* Measure potential on-line? */
	   int HQPotItr,		/* Iterations per potential measurement */
	   Real FacSmea,                  /* Smearing factor */
	   int NumSmea,		/* Number of smearing iterations */

	   int GFixP,			/* Fix to Landau gauge? */
	   int GFixItr,		/* Iterations per Landau gauge fixing */
	   int ORlxDo,		/* Do Over-relaxation in gauge fixing */
	   Real OrPara,			/* Over-relaxation parameter in gauge-fixing */
	   )
{
  /* Include any hooks in the start macro */
  START_CODE();			/*# May include declarations and/or code. */

  multi1d<LatticeColorMatrix> u_tmp(Nd); /* temp. gauge field cfg. */
  multi1d<LatticeColorMatrix> u_sm(Nd); /* smeared gauge field */

  multi1d<DComplex> pollp(Nd);	/* Polyakov loop */
  Double w_plaq;		/* Whole plaquette */
  Double s_plaq;		/* Spatial plaquette */
  Double t_plaq;		/* Temporal (thermal) plaquette */
  Double link;			/* Trace of link (not gauge invariant) */

  int GFMax;		/* Maximum number of gauge-fixing relaxations */
  Real GFAccu;			/* Gauge-fixing relaxation accuracy */
  int nrl_gf;

  /* Start executable code */
  NumTry = 0;
  NumFail = 0;

  /* hardwire some gauge fixing parameters */
  GFAccu = 1.0e()-6;
  GFMax = 1000;

  /* Perform Monte Carlo algorithm and measurements, as requested */

  /* Do NumItr iterations of NOver overrelaxation steps followed by */
  /* one heatbath step */
  for(int n_loop = 1; n_loop <= NumItr; ++n_loop)
  {
    TotalItr = TotalItr + 1;
    int ItrNumber = TotalItr;         /* Want a different name here. */
    
    /* This is the very first namelist group printed in an iteration */
    push(xml_out,"NewIteration");
    write(xml_out, "ItrNumber", ItrNumber);
    pop(xml_out);

    /* Perform a Monte Carlo iteration */
#if 1
#warning "NEED MORE GENERAL SET USAGE HERE"
    mciter(u, NOver, NHeat, NumTry, NumFail, rb, neighsubl);
#endif
    
    /* Carry out measurements every MesItr iterations */
    if ((TotalItr % MesItr) == 0 )
    {
      MesPlq(u, w_plaq, s_plaq, t_plaq, link);

      /* Measure Polyakov loops */
      for(int mu = 0; mu < Nd; ++mu)
	polylp(u, pollp[mu], mu);

      push(xml_out,"obsvbl1");
      write(xml_out, "TotalItr", TotalItr);
      write(xml_out, "w_plaq", w_plaq);
      write(xml_out, "s_plaq", s_plaq);
      write(xml_out, "t_plaq", t_plaq);
      write(xml_out, "link", link);
      write(xml_out, "pollp", pollp);
      pop(xml_out);
      push(xml_out,"obsvbl2");
      write(xml_out, "NumTry", NumTry);
      write(xml_out, "NumFail", NumFail);
      pop(xml_out);

      if ( GlueP == YES )
	fuzglue (u, t_dir, BlkAccu, BlkMax);

      if ( WilsonP == YES ) 
	wilslp (u, t_dir, 3);

      /* Print statistics every PrtMes measurements */
      if ((TotalItr & MesItr*PrtMes) == 0)
      {
	FPRINTF(trm_out,"\n     Heatbath acceptance =%14.7g\n",
		TO_REAL(NumTry-NumFail)/TO_REAL(NumTry));

      } /* end of on-line statistics display */
    } /* end of measurements */
    
    /* Calculate topological charge and the action/continuum instanton action ratio. */
    if ((TotalItr & TopolItr) == 0)
    {
      /* Topology */
      if( TopolP == YES  )
      {
	u_tmp = u;
	topol(u_tmp, TopAccu, ActAccu, NumCool, NumTop);
      }
    }
    
    /* Calculate topological charge and the action/continuum instanton action ratio. */
    if ((TotalItr & HQPotItr) == 0)
    {
      if ( HQPotP )   /* First regular potl on iso or aniso lattices: */
      {
	u_sm = u;

	u_tmp[t_dir] = u_sm[t_dir];

	/* Smear the space-like links NumSmea times */
	for(int i = 1; i <= NumSmea; ++i)
	{
	  for(int mu = 0; mu  < Nd; ++mu)
	    if( mu != t_dir )
	      smear(u_sm, u_tmp[mu], mu, 0, FacSmea, BlkAccu, BlkMax, t_dir);

	  u_sm = u_tmp;
	}

	
	push(xml_out,"Smeared_Wilson_Loops");
	write(xml_out, "NumSmea", NumSmea);
	write(xml_out, "FacSmea", FacSmea);
	write(xml_out, "t_dir", t_dir);
	pop(xml_out);
	wilslp (u_sm, t_dir, 6);


	/* In ANISO case do sideways potl ONLY if j_decay != t_dir : */
  
	if( AnisoP == YES && j_decay != t_dir && j_decay < Nd ) 
	{
	  u_sm = u;
	  u_tmp[j_decay] = u_sm[j_decay];

	  /* Smear the space-like links NumSmea times */
	  for(int i = 1; i <= NumSmea; ++i)
	  {
	    for(int mu = 0; mu  < Nd; ++mu)
	      if( mu != j_decay )
		smear(u_sm, u_tmp[mu], mu, 0, FacSmea, BlkAccu, BlkMax, j_decay);
	  }
	  u_sm = u_tmp;

	  push(xml_out,"Smeared_Wilson_Loops");
	  write(xml_out, "NumSmea", NumSmea);
	  write(xml_out, "FacSmea", FacSmea);
	  write(xml_out, "j_decay", j_decay);
	  pop(xml_out);
	  wilslp(u_sm, j_decay, 6);
	}               /*  end of sideways potl  */

      }                 /*  closes HQPotP == YES  */
    }
    
    /* Fix to Landau gauge to measure Landau gauge link */
    if ((TotalItr & GFixItr) == 0)
    {
      if( GFixP == YES  )
      {
	u_tmp = u;
	gfix(u_tmp, Nd, GFAccu, GFMax, nrl_gf, ORlxDo, OrPara);
      }
    }
    
    /* This is absolutely the last namelist group Printed in this trajectory */
    push(xml_out,"EndIteration");
    write(xml_out, "ItrNumber", ItrNumber);
    pop(xml_out);
  }                     /* end loop over iterations */
  
  /* Close out any other code */
  END_CODE();
}

}  // end namespace Chroma
