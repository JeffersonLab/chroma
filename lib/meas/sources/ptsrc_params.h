// -*- C++ -*-
// $Id: ptsrc_params.h,v 1.1 2005-10-28 21:06:41 edwards Exp $
/*! \file
 *  \brief Point source parameters
 *
 * This code is common to different fermions and propagators
 */

#ifndef __ptsrc_params_h__
#define __ptsrc_params_h__

namespace Chroma
{

  //! Point source parameters
  /*! @ingroup sources */
  struct PointSourceConstParams
  {
    PointSourceConstParams() {}
    PointSourceConstParams(XMLReader& in, const std::string& path);
    
    multi1d<int>     t_source;        /*!< source location */
  };


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, PointSourceConstParams& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const PointSourceConstParams& param);

}  // end namespace Chroma


#endif
