// -*- C++ -*-
// $Id: monte.h,v 3.0 2006-04-03 04:59:07 edwards Exp $
/*! \file
 *  \brief Routine for doing the hybrid (monte carlo) algorithm. 
 */

#ifndef __monte_h__
#define __monte_h__

namespace Chroma {

//! Routine for doing the hybrid (monte carlo) algorithm.
/*!
 * \ingroup heatbath
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
	   );

}  // end namespace Chroma

#endif
