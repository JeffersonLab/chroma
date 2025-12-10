/* $Id: packer_p4_pad.c,v 1.1 2007-09-12 19:33:13 bjoo Exp $ */

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

#if SSE_PRECISION == 32
#warning "Building Packer for 32 Bits"
typedef float u_mat_array[3][3][2];
#elif SSE_PRECISION == 64
#warning "Building Packer for 64 Bits"
typedef double u_mat_array[3][3][2];
#else
#error "Precision Not supported, define SSE_PRECISION"
#endif

/* here's how packing works  ...pack_gague_field...for the 32-bit parallel implemenation of the dslash
   u_tmp <-- u
   in lexical order of u(site, mu) the packed gauge fields are laid out as follows:
   u_tmp(site,0)
   u_tmp(site+1,0)
   u_tmp(site,1)
   u_tmp(site+1,1)
   u_tmp(site,2)
   u_tmp(site+1,2)
   u_tmp(site,3)
   u_tmp(site+1,3)
*/
 
 
/*no_funnystuff_pack_gauge_field does not mix directions and sites as pack_gague_field does, and must be used
for all of the other P4 dslashes */


void pack_gauge_field(int volume, u_mat_array *u, u_mat_array *u_tmp) /* pass &u[0][0][0] */
{
  int ix,mu,cb;
  u_mat_array v[8];

  for (cb=0;cb<2;cb++)
  {
    for (ix=0;ix<volume;ix+=2)
    {
      for (mu=0;mu<4;mu++)
      {
	memcpy(v + 2*mu, u+ix + volume*(cb + 2*(mu)), sizeof(u_mat_array)); /*gauge_field[ix][mu]*/
	memcpy(v + 2*mu+1, u+(ix+1) + volume*(cb + 2*mu), sizeof(u_mat_array)); /*gauge_field[ix+1][mu]*/
      }		

      for (mu=0;mu<4;mu++)
      {
	memcpy(u_tmp+mu+4*(ix+volume*(cb)), v + mu, sizeof(u_mat_array));
	memcpy(u_tmp+mu+4*(ix+1+volume*(cb)), v + 4+mu, sizeof(u_mat_array));
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
	memcpy(v + mu, u_tmp+mu+4*(ix + volume*(cb)), sizeof(u_mat_array));
	memcpy(v + 4+mu,u_tmp+mu+4*(ix+1+ volume*(cb)) , sizeof(u_mat_array));
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
