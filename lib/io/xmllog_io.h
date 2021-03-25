// -*- C++ -*-
/*! \file
 * \brief Singleton instances of xml output
 */

#ifndef _XMLLOG_IO_H_
#define _XMLLOG_IO_H_

#include "chromabase.h"
#include "singleton.h"


namespace Chroma 
{

  //! XML output holder
  /*! \ingroup io */
  typedef Chroma::SingletonHolder<XMLFileWriter> TheXMLOutputWriter;

  //! XML log holder
  /*! \ingroup io */
  typedef Chroma::SingletonHolder<XMLFileWriter> TheXMLLogWriter;

} // End namespace Chroma

#endif
