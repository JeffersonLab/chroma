// $Id: common_declarations.h,v 1.4 2003-03-29 05:21:48 edwards Exp $

#ifndef COMMON_DECLS_INCLUDE
#define COMMON_DECLS_INCLUDE

#if defined(MAIN)
#define EXTERN
#else
#define EXTERN extern
#endif

/* More variables describing state of system and algorithm */
EXTERN int FermTypeP; /* Type of Fermions use, Wilson or Staggered */
EXTERN int MaxCG;
EXTERN int FermiP;
EXTERN Real Nf;          /* Number of Flavours      */
EXTERN Real BetaMC;
EXTERN Real BetaMD;
EXTERN Real KappaMC;   /* Input for Wilson, or = 1/2MassMC for Staggered */
EXTERN Real KappaMD;   /* Input for Wilson, or = 1/2MassMD for Staggered */
EXTERN Real MassMC;    /* Used only for Staggered fermions */
EXTERN Real MassMD;    /* Used only for Staggered fermions */
EXTERN Real dt;
EXTERN Real RsdCGMC;
EXTERN Real RsdCGMD;
EXTERN Real RsdCGMDMax;
EXTERN Real RsdCGMDFctr;
EXTERN int IntrplOrd;  /* Order of interpolation */
EXTERN int Algorithm;  /* No Metropolis step (i.e., use Hybrid */
EXTERN int AlgETrj;    /* Specify trajectory length algorithm */
EXTERN int AlgLPStp;   /* Algorithm requires P and psi at end of traj. */
EXTERN int RefNextTrj; /* Refresh P or chi for next trajectory? */
EXTERN int RalgP;      /* R-algorithm mode? */
EXTERN Real LamPl;         /* For stochastic accept/reject test */
EXTERN Real LamMi;         /* For stochastic accept/reject test */
EXTERN Real AlpLog;        /* For stochastic estimate of Log(1+x) */
EXTERN Real AlpExp;        /* For stochastic estimate of Exp(x)   */
EXTERN int MonitordH;  /* Measure dH for every step on the trajectory */
EXTERN Real tau0;          /* Average trajectory length */

/* New variables related to Action improvement */
extern int InvType;  /* Type of fermion inverter */
extern int GaugeAct; /* Type of Gauge action (Wilson,Symanzik,...) */
extern int GlueImp;  /* Level of Symanzik improvement in dsdu */
extern Real GlueCoeffRT;   /* Coefficient of 1x2 rectangle in dsdu */
extern Real GlueCoeffPG;   /* Coefficient of parallelogram in dsdu */
extern int FermAct; /* Type of Fermion action (Wilson,Clover,...) */
extern Real MRover;   /* Clover inverter over-relaxtion parameter */
extern Real ClovCoeff;       /* Relative Clover coefficient (=1 at tree level) */
extern Real u0;         /* Tadpole improvement factor  (Tr U_p/3)^(1/4) */
extern Real H_parity;        /* Strenght of parity breaking term */
extern Real ClovCoeffR; /* Spatial Clover coefficient */
extern Real ClovCoeffT; /* Temporal Clover coefficient */


/* variable which controls number of namelist sites */
// EXTERN unsigned global_namelist_vol;

// Initialization section
#if defined(MAIN)
int FftInitP =0;
int SchrFun=0;
int SftInitP=0;

int InvType=0;  /* Type of fermion inverter */
int GaugeAct=0; /* Type of Gauge action (Wilson,Symanzik,...) */
int GlueImp=0;  /* Level of Symanzik improvement in dsdu */
int FermAct=0; /* Type of Fermion action (Wilson,Clover,...) */
Real MRover=1.0;   /* Clover inverter over-relaxtion parameter */
Real u0=0.0;         /* Tadpole improvement factor  (Tr U_p/3)^(1/4) */
Real ClovCoeff=0.0;       /* Relative Clover coefficient (=1 at tree level) */
Real ClovCoeffR=0.0; /* Spatial Clover coefficient */
Real ClovCoeffT=0.0; /* Temporal Clover coefficient */
Real H_parity;        /* Strenght of parity breaking term */
int PolyDegHalf=0;  /* degree(p)/2                         */
int RatPolyDeg=0;  /* Degree of polynomials in ratio approx  */

int OverAuxAct=UNPRECONDITIONED_WILSON;   /* Kernel of overlap,dwf */
Real OverMass=0.0;   /* Overlap mass */
int NOperEig=0;       /* Number of eigenvectors kept */
int NOperEigDim=0;    /* Dimension of eigenvectors kept */
#endif

#endif
