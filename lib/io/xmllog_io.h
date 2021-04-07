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
  struct TheXMLOutputWriterInstance {};
  typedef Chroma::SingletonHolder<XMLFileWriter, TheXMLOutputWriterInstance> TheXMLOutputWriter;

  //! XML log holder
  /*! \ingroup io */
  struct TheXMLLogWriterInstance {};
  typedef Chroma::SingletonHolder<XMLFileWriter, TheXMLLogWriterInstance> TheXMLLogWriter;

} // End namespace Chroma

#endif
