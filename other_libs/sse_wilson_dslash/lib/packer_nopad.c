/* $Id: packer_nopad.c,v 1.1 2007-09-12 19:33:13 bjoo Exp $ */

/*  does the packing of the gague fields needed by the dslash */

/* this must be called before the P4 dslash is called */
/* pack_gauge_field(int volume, u_mat_array *u, u_mat_array *u_tmp) */
/* volume: # of sites */
/* u: input normal SZIN ordered gauge field */
/* u_tmp: output packed gague field for use by P4 dslash */

#include <sse_config.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef DMALLOC
#include <stdlib.h>
#else
#include <dmalloc.h>
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>

typedef SSEREAL u_mat_array[3][3][2];

void pack_gauge_field(int volume, u_mat_array *u, u_mat_array *u_tmp) /* pass &u[0][0][0] */
{
  int ix,mu,cb;
  for (cb=0;cb<2;cb++)
  {

    for (ix=0;ix<volume;ix++)
    {
      for (mu=0;mu<4;mu++)
      {
	memcpy(u_tmp+ (mu + 4*(ix+volume*(cb))), 
	       u+ix + volume*(cb + 2*(mu)), 
	       sizeof(u_mat_array)); /*gauge_field[ix][mu]*/
      }		
    }
  }
}

void unpack_gauge_field(int volume, u_mat_array *u_tmp, u_mat_array *u)  /* pass &u[0][0][0] */
{
  int ix,mu,cb;
  u_mat_array v[8];
  for (cb = 0; cb<2; cb++)
  {
    for (ix=0;ix<volume;ix+=2)
    {
      for (mu=0;mu<4;mu++)
      {
	memcpy(v + 2*mu, u_tmp+ix + volume*(cb + 2*mu), sizeof(u_mat_array));
	memcpy(v + 2*mu+1, u_tmp+(ix+1) + volume*(cb + 2*mu), sizeof(u_mat_array));
      }

      for (mu=0;mu<4;mu++)
      {
	memcpy(u+ix + volume*(cb + 2*mu), v + 2*mu, sizeof(u_mat_array));
	memcpy(u+(ix+1) + volume*(cb + 2*mu), v + 2*mu+1, sizeof(u_mat_array));
      }
    }
 
  }
}

#ifdef __cplusplus
}
#endif
