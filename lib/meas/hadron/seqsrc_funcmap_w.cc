// $Id: seqsrc_funcmap_w.cc,v 2.1 2005-09-26 04:48:35 edwards Exp $
/*! \file
 *  \brief Sequential source function map
 */

#include "meas/hadron/seqsrc_funcmap_w.h"
#include "meas/hadron/barseqsrc_w.h"
#include "meas/hadron/mesonseqsrc_w.h"

namespace Chroma
{
 
  namespace SeqSourceCallMapEnv
  { 
    bool registerAll(void) 
    {
      bool success = true;
      success &= MesonSeqSourceCallMapEnv::registered;
      success &= BaryonSeqSourceCallMapEnv::registered;
      return success;
    }

    bool registered = registerAll();
  }

}
