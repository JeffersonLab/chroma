// $Id: common_io.h,v 1.2 2003-04-17 20:10:30 dgr Exp $

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
  int sink_smearingparam;
  int sink_type;
  int sink_direction;
};



#endif
