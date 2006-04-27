// $Id: inline_io_aggregate.cc,v 3.1 2006-04-27 02:35:38 edwards Exp $
/*! \file
 *  \brief Inline IO aggregator
 */

#include "meas/inline/io/inline_io_aggregate.h"
#include "meas/inline/io/inline_qio_write_obj.h"
#include "meas/inline/io/inline_qio_write_erase_obj.h"
#include "meas/inline/io/inline_qio_read_obj.h"
#include "meas/inline/io/inline_erase_obj.h"
#include "meas/inline/io/inline_list_obj.h"
#include "meas/inline/io/inline_szin_read_obj.h"
#include "meas/inline/io/inline_szin_write_obj.h"
#include "meas/inline/io/inline_nersc_read_obj.h"
#include "meas/inline/io/inline_nersc_write_obj.h"

#include "meas/inline/io/inline_xml_write_obj.h"

#include "meas/inline/io/inline_gaussian_obj.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineIOAggregateEnv
  {
    bool registerAll() 
    {
      bool success = true; 

      // Tasks
      success &= InlineQIOReadNamedObjEnv::registered;
      success &= InlineQIOWriteNamedObjEnv::registered;
      success &= InlineQIOWriteEraseNamedObjEnv::registered;
      success &= InlineEraseNamedObjEnv::registered;
      success &= InlineListNamedObjEnv::registered;

      success &= InlineGaussianInitNamedObjEnv::registered;

      success &= InlineSZINReadNamedObjEnv::registered;
      success &= InlineSZINWriteNamedObjEnv::registered;

      success &= InlineNERSCReadNamedObjEnv::registered;
      success &= InlineNERSCWriteNamedObjEnv::registered;

      success &= InlineXMLWriteNamedObjEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
