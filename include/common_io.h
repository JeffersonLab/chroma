// $Id: common_io.h,v 1.3 2003-08-28 22:06:32 ikuro Exp $

#ifndef __common_io_h__
#define __common_io_h__

#if defined(MAIN)
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN BinaryReader cfg_in;
EXTERN BinaryWriter cfg_out;
// EXTERN NmlReader nml_in;
EXTERN NmlWriter nml_out;
// EXTERN FILE * trm_in;
// EXTERN FILE * trm_out;
// EXTERN int Binary_io_location;

//  Finally, some (temporary) routines associated with the headers

struct PropHead{
  Real kappa;
  int source_smearingparam;
  int source_type;		// S-wave (0) or P-wave  (1)
  int source_direction;         // S-wave (0);   P-wave x(0) y(1) z(2)
  int source_laplace_power;
  int sink_smearingparam;
  int sink_type;
  int sink_direction;
  int sink_laplace_power;
};



#endif
