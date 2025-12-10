/*******************************************************************************
 * $Id: sse_align.h,v 1.2 2007-09-12 20:29:39 bjoo Exp $
 * 
 *
 * Define alignment macro ALIGN. This can vary depending on compiler 
 *
 * Author: Balint Joo 
 * Date: 11/13/2002
 *
 *******************************************************************************/

/* Include guard... */
#ifndef __INCLUDED_SSE_ALIGN_H__
#define __INCLUDED_SSE_ALIGN_H__

#include <sse_config.h>

#ifndef ALIGN

/* Gnu compilers */
#if __GNUC__

/* We no longer support GCC v2 */
#if __GNUC__ == 2
#error "GNU C Version 2 no longer supported"
#endif 

#if __GNUC__ == 3 && __GNUC_MINOR__ <= 3
#error "If Using GCC v3, use version GCC v3.3 or better"
#endif
 
#define ALIGN __attribute__ ((aligned (16)))

#else  /* __GNUC__ */

/* Define it as empty for unknown compiler */
#define ALIGN 

#endif /* if __GNUC__ */


#endif /* ALIGN */



/* End of include guard */
#endif 


