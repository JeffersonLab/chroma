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

#include "meas/inline/io/inline_rng.h"


#include "meas/inline/io/inline_xml_write_obj.h"

#include "meas/inline/io/inline_gaussian_obj.h"
#include "meas/inline/io/inline_usqcd_read_ddpairs_prop.h"
#include "meas/inline/io/inline_usqcd_write_ddpairs_prop.h"
#include "meas/inline/io/inline_milc_write_stag_source.h"

#include "meas/inline/io/inline_eigen_bin_colvec_read_obj.h"
#include "meas/inline/io/inline_eigen_lime_colvec_read_obj.h"
#include "meas/inline/io/inline_eigen_bin_lime_colvec_read_obj.h"
#include "meas/inline/io/inline_read_map_obj_disk.h"
#include "meas/inline/io/inline_copy_map_obj.h"
#include "meas/inline/io/inline_write_timeslice_map_obj_disk.h"

#include "chroma_config.h"
#ifdef BUILD_QOP_MG
#include "meas/inline/io/inline_erase_mg_space.h"
#endif

#ifdef BUILD_QUDA
#include "meas/inline/io/inline_erase_quda_multigrid_space.h"
#endif

#ifdef BUILD_MGPROTO
#include "meas/inline/io/inline_erase_mg_proto_space.h"
#ifdef BUILD_QPHIX
#include "meas/inline/io/inline_erase_mg_proto_qphix_space.h"
#endif
#endif

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
	success &= InlineSetRNGEnv::registerAll();

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

	success &= InlineEigenBinColVecReadNamedObjEnv::registerAll();
	success &= InlineEigenLimeColVecReadNamedObjEnv::registerAll();
	success &= InlineEigenBinLimeColVecReadNamedObjEnv::registerAll();

	success &= InlineXMLWriteNamedObjEnv::registerAll();


	// QIO USQCD DD PAIRS Reader
	success &= InlineUSQCDReadDDPairsPropEnv::registerAll();
	success &= InlineUSQCDWriteDDPairsPropEnv::registerAll();

	// MILC Source Writer
	success &= InlineMILCWriteStagSourceEnv::registerAll();

	// MapObjDisk reader
       	success &= InlineReadMapObjDiskEnv::registerAll();
       	success &= InlineCopyMapObjEnv::registerAll();
       	success &= InlineWriteTimeSliceMapObjDiskEnv::registerAll();

#ifdef BUILD_QOP_MG
	success &= InlineEraseMGSpaceEnv::registerAll();
#endif
	
#ifdef BUILD_QUDA
	success &= InlineEraseQUDAMULTIGRIDSpaceEnv::registerAll();
#endif

#ifdef BUILD_MGPROTO
	success &= InlineEraseMGProtoSpaceEnv::registerAll();
#ifdef BUILD_QPHIX
	success &= InlineEraseMGProtoQPhiXSpaceEnv::registerAll();
#endif
#endif

	registered = true;
      }
      return success;
    }
  }

}
