// $Id: source_const_aggregate.cc,v 3.8 2009-06-24 19:55:40 jbulava Exp $
/*! \file
 *  \brief All make source constructors
 */

#include "meas/sources/source_const_aggregate.h"

#include "meas/sources/pt_source_const.h"
#include "meas/sources/sh_source_const.h"
#include "meas/sources/norm_sh_source_const.h"
#include "meas/sources/wall_source_const.h"
#include "meas/sources/mom_source_const.h"
#include "meas/sources/partwall_source_const.h"

#include "meas/sources/rndz2wall_source_const.h"
#include "meas/sources/rndzNwall_source_const.h"
#include "meas/sources/dilutezN_source_const.h"
#include "meas/sources/dilute_zN_eigvec_source_const.h"
#include "meas/sources/diluteGrid_source_const.h"

#include "meas/sources/sf_pt_source_const.h"
#include "meas/sources/sf_sh_source_const.h"
#include "meas/sources/sf_wall_source_const.h"

namespace Chroma
{

  //! Registration aggregator
  namespace QuarkSourceConstructionEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Sources
	success &= PointQuarkSourceConstEnv::registerAll();
	success &= ShellQuarkSourceConstEnv::registerAll();
	success &= NormShellQuarkSourceConstEnv::registerAll();
	success &= WallQuarkSourceConstEnv::registerAll();
	success &= RandZ2WallQuarkSourceConstEnv::registerAll();
	success &= RandZNWallQuarkSourceConstEnv::registerAll();
	success &= MomWallQuarkSourceConstEnv::registerAll();
	success &= PartialWallQuarkSourceConstEnv::registerAll();
	success &= DiluteZNQuarkSourceConstEnv::registerAll();
	success &= DiluteZNEigVecQuarkSourceConstEnv::registerAll();

	success &= DiluteGridQuarkSourceConstEnv::registerAll();

	success &= SFPointQuarkSourceConstEnv::registerAll();
	success &= SFShellQuarkSourceConstEnv::registerAll();
	success &= SFWallQuarkSourceConstEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
