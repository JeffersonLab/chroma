// $Id: inline_io_aggregate.cc,v 3.3 2008-04-25 18:09:01 bjoo Exp $
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
#include "meas/inline/io/inline_usqcd_read_ddpairs_prop.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineIOAggregateEnv
  {
    namespace
    {
      //! Local registration flag
      bool registered = false;
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Tasks
	success &= InlineQIOReadNamedObjEnv::registerAll();
	success &= InlineQIOWriteNamedObjEnv::registerAll();
	success &= InlineQIOWriteEraseNamedObjEnv::registerAll();
	success &= InlineEraseNamedObjEnv::registerAll();
	success &= InlineListNamedObjEnv::registerAll();

	success &= InlineGaussianInitNamedObjEnv::registerAll();

	success &= InlineSZINReadNamedObjEnv::registerAll();
	success &= InlineSZINWriteNamedObjEnv::registerAll();

	success &= InlineNERSCReadNamedObjEnv::registerAll();
	success &= InlineNERSCWriteNamedObjEnv::registerAll();

	success &= InlineXMLWriteNamedObjEnv::registerAll();


	// QIO USQCD DD PAIRS Reader
	success &= InlineUSQCDReadDDPairsPropEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
