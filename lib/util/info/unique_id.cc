// $Id: unique_id.cc,v 1.2 2007-02-25 22:37:35 edwards Exp $
/*! \file
 *  \brief Generate a unique id
 */

#include "chromabase.h"
#include "util/info/unique_id.h"
#include <time.h>

namespace Chroma 
{
  //! Return a unique id
  /*!
   * \ingroup info
   *
   *  The id is return type a string. This gives maximal flexibility allowing
   *  the way the ID is generated to change in the future.
   */
  std::string uniqueId()
  {
    START_CODE();

    // The id is the seconds since RGE started work at JLab...
    // Chroma certainly did not live before this date.
    struct tm Oct_1_struct = {0};
    time_t Oct_1_t;
    Oct_1_struct.tm_year = 99;
    Oct_1_struct.tm_mon  = 10;
    Oct_1_struct.tm_mday = 1;
    Oct_1_t = mktime(&Oct_1_struct);
    if (Oct_1_t == (time_t)-1)
    {
      QDPIO::cerr << __func__ << ": some error generating ID" << endl;
      exit(1);
    }

    double dd = difftime(time(NULL), Oct_1_t);
    unsigned long int ld = (unsigned long int)(dd);
    std::ostringstream foo;
    foo << ld;

    END_CODE();

    return foo.str();
  }

}  // end namespace Chroma
