// $Id: primitives.h,v 2.0 2005-09-25 21:04:25 edwards Exp $

#error "OBSOLETE - DO NOT USE. ONLY FOR REFERENCE"

#ifndef PRIMITIVES_INCLUDE
#define PRIMITIVES_INCLUDE

/*
 * Symbolic names
 */
#define F_SWAP  0		 /* DUPLEX indices */
#define B_SWAP  1

/*
 * Valid option values
 */
#define PURE_GAUGE  0
#define WILSON_GAUGE  0
#define SYMANZIK_GAUGE  1
#define MANTON_GAUGE  2
#define WILSON_FERMIONS  1
#define STAGGERED_FERMIONS  2
#define WILSON  0
#define PARITY_BREAKING_WILSON  1
#define CLOVER  2
#define PERFECT_WILSON_FERMIONS  3
#define W12_FERMIONS  4
#define UNPRECONDITIONED_WILSON  5
#define UNPRECONDITIONED_CLOVER  6
#define OVERLAP_POLE  7
#define OVERLAP_INVERSE  8
#define TRUNC_OVERLAP  9
#define ZOLOTAREV_4D  12
#define ZOLOTAREV_5D  14
#define OVERLAP_5D  15
#define OVERLAP_DWF  16
#define DWF  17
#define OVERLAP_DWF_4D  18
#define DWF_4D  19
#define DWF_TRANSF  20
#define DWF_POLE  21
#define TRUNC_DWF_POLE  22
#define PRECONDITIONED_DWF  23
#define EXTENDED_OVERLAP  24
#define PRECONDITIONED_EXTENDED_OVERLAP  25
#define SMEARED_LAPLACIAN_WILSON  26
#define PLANAR_WILSON  27
#define HAMBER_WU  28
#define PROJECTED_PRECONDITIONED_DWF  29

#define WILSON_DSLASH  501
#define DWF_DSLASH  502

#define STAGGERED  10
#define NAIK  11

#define QUADRATIC_BOSON  30

// #define CG_INVERTER  21
// #define MR_INVERTER  22
// #define BICG_INVERTER  23
// #define CR_INVERTER  24
#define HMD  91
#define HMC  92
#define HMDC  94
#define HMCC  95
#define PHMD  96
#define PHMC  98
#define PHMCN  99
#define RHMD  996
#define RHMC  998
#define RHMCN  999
#define KPHB  201
#define CrHB  202
#define FIXED_LENGTH  111
#define EXPONENTIAL_LENGTH  113
#define SCHROEDINGER_BACKGROUND  2
#define HOT  1
#define COLD  0
#define RUN_CONT  -1
#define STAT_CONT  2
#define STAT_RESTART  3
#define FORWARD  1
#define BACKWARD  -1
#define YES  1
#define NO  0
#define BE_BINARY_LOCATION  0
#define FE_BINARY_LOCATION  1
#define OPTION_REPLACE  11
#define OPTION_NEGATE  13
#define OPTION_ADD  17
#define OPTION_SUBTRACT  19
#define OPTION_REAL_PART  23
#define OPTION_IMAGINARY_PART  29
#define OPTION_COMPLEX_PART  31
#define OPTION_REUNITARIZE  37
#define OPTION_REUNITARIZE_ERROR  41
#define OPTION_REUNITARIZE_LABEL  43
#define OPTION_TWELTH_ORDER  47
#define OPTION_EXACT  53
#define OPTION_POINT_SOURCE  59
#define OPTION_WALL_SOURCE  61
#define OPTION_POINT_SINK  67
#define OPTION_WALL_SINK  71
#define OPTION_POINT_AND_WALL_SINK  73
#define OPTION_SHELL_SOURCE  79
#define OPTION_BNDST_SOURCE  83
#define OPTION_POINT_AND_BNDST_SOURCE  89
#define OPTION_SHELL_AND_BNDST_SOURCE  97
#define OPTION_POINT_AND_SHELL_AND_BNDST_SOURCE  101
#define OPTION_SHELL_SINK  103
#define OPTION_POINT_AND_SHELL_SINK  107
#define OPTION_BNDST_SINK  113
#define OPTION_POINT_AND_BNDST_SINK  127
#define OPTION_SHELL_AND_BNDST_SINK  131
#define OPTION_POINT_AND_SHELL_AND_BNDST_SINK  137
#define OPTION_WALL_WVF  191
#define OPTION_DELTA_WVF  179
#define OPTION_PWV_DELTA_WVF  181
#define OPTION_DWV_DELTA_WVF  187
#define OPTION_GAUGE_INV_GAUSSIAN_WVF  189
#define OPTION_PWV_GAUGE_INV_GAUSS_WVF  193
#define OPTION_DWV_GAUGE_INV_GAUSS_WVF  197
#define OPTION_GAUSSIAN_WVF  139
#define OPTION_PWV_GAUSSIAN_WVF  151
#define OPTION_DWV_GAUSSIAN_WVF  163
#define OPTION_EXPONENTIAL_WVF  149
#define OPTION_PWV_EXPONENTIAL_WVF  157
#define OPTION_DWV_EXPONENTIAL_WVF  167
#define OPTION_WUPPERTAL_WVF  199
#define OPTION_PWV_WUPPERTAL_WVF  203
#define OPTION_DWV_WUPPERTAL_WVF  209

#endif
