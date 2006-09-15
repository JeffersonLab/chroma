// -*- C++ -*-
// $Id: xmllog_io.h,v 3.1 2006-09-15 02:48:21 edwards Exp $
/*! \file
 * \brief Singleton instances of xml output
 */

#ifndef _XMLLOG_IO_H_
#define _XMLLOG_IO_H_

#include "chromabase.h"
#include "singleton.h"


namespace Chroma 
{

  /* typedef SingletonHolder< XMLReader > TheXMLInputReader; */

  // Create each singleton with a slightly different lifetime policy
  // This trick is used to disambiguate only

  //! XML output holder
  /*! \ingroup io */
  typedef SingletonHolder<XMLFileWriter, CreateUsingNew,
			  DefaultLifetime1,
			  SingleThreaded> TheXMLOutputWriter;

  //! XML log holder
  /*! \ingroup io */
  typedef SingletonHolder<XMLFileWriter, CreateUsingNew,
			  DefaultLifetime2,
			  SingleThreaded> TheXMLLogWriter;

} // End namespace Chroma

#endif
