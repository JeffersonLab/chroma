/* In the original version of this file, the QDP precision macro was defined in the compiler command line.
 * I don't know how to make the 2 objects from the same source file with automake. So I am wrapping them
 * up here, and setting the precision macro. The nice part is that I don't need separate compile commands
 */

#undef QDP_Precision
#define QDP_Precision 'D'
/* All the work is in here */
#include "./wilsonmg_p.c"
#undef QDP_Precision

