#ifndef SHIFT_TABLES_PARSCALAR_H
#define SHIFT_TABLES_PARSCALAR_H

#include <stddef.h>

#include <sse_config.h>

#if SSE_PRECISION == 32
#include <types32.h>
#else
#include <types64.h>
#endif 

#ifdef __cplusplus
extern "C" {
#endif

  typedef enum {
    DECOMP_SCATTER=0,
    DECOMP_HVV_SCATTER,
    RECONS_MVV_GATHER,
    RECONS_GATHER
  } HalfSpinorOffsetType;


  // int getSubgridVol();
  int getSubgridVolCB();
  int getSubgridVol3D();
  void make_shift_tables(int bound[2][4][4], 
			 halfspinor_array* chi1, 
			 halfspinor_array* chi2,
			 halfspinor_array* recv_bufs[2][4],
			 halfspinor_array* send_bufs[2][4],
			 void (*getSiteCoords)(int coord[], int node, int linearsite), 
			 int (*getLinearSiteIndex)(const int coord[]),
			 int (*getNodeNumber)(const int coord[]));

  void make_shift_tables_3d(int bound[2][4][3],
			    void (*getSiteCoords)(int coord[], int node, int linearsite), 
			    int (*getLinearSiteIndex)(const int coord[]),

			    int (*getNodeNumber)(const int coord[]));


  void free_shift_tables(void);
  void free_shift_tables_3d(void);

  
  /*! This is the key routine. It is a table lookup for an offset 
    into the half spinor array. It indirects us to the tail regions
    as appropriate for sites to be communicated off-node 
    
    \param type (HalfSpinorOffsetType). The type of operation whose offset we want (eg: DECOMP_SCATTER, DECOMP_HVV_SCATTER, RECONS_MVV_GATHER, RECONS_GATHER)
    \param site (int) The site we want to do the offset lookup for
    \param mu   (int) The index for the gamma matrix in the projector (0-3)
  */
  /* int halfspinor_buffer_offset(HalfSpinorOffsetType type, int site, int mu) ;  */

#ifdef __cplusplus
};
#endif


#endif
