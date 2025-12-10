#ifndef SSE_DSLASH_H
#define SSE_DSLASH_H

#include "sse_config.h"

#ifdef __cplusplus
  extern "C" {
#endif
  
    /* The C Function declarations */
    /*! Initialize the SU3 Dslash */
    /*! \param latt_size  An integer array speficying the size of the lattice */
    void init_sse_su3dslash(const int* latt_size,
			    void (*getSiteCoords)(int coord[], int node, int linearsite),
			    
			    int (*getLinearSiteIndex)(const int coord[]),
			    int (*nodeNumber)(const int coord[]));

    /*! Finalize the SU3 Dslash */
    void free_sse_su3dslash(void);

    /*! Apply the SU3 Dslash */
    /*! \param u  Pointer to the packed gauge field. 
        \param psi Pointer to the source vectors first element (regardless of cb)
	\param chi Pointer to the result vectors first element (regardless of cb)
	\param isign -1 for applying Dslash Dagger +1 for applying Dslash 
	
	\param cb  The checkerboard index of the source std::vector.

        NB: psi and chi extend over the whole lattice in this apparoach (they are not solely one checkerboards worth of std::vector but a full latticeworth. Offsets are worked out internally.
    */
    void sse_su3dslash_wilson(SSEREAL* u, SSEREAL *psi, SSEREAL *res, int isign, int cb);
  

    /*! Pre post receives for bigger ops? */
    void sse_su3dslash_prepost_receives(void);

#ifdef __cplusplus
  };
#endif

#endif
