// $Id: common_io.h,v 1.1 2003-04-07 04:48:19 edwards Exp $

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

#endif
