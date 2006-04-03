// -*- C++ -*-
// $Id: xmllog_io.h,v 3.0 2006-04-03 04:58:56 edwards Exp $
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

  typedef SingletonHolder< XMLFileWriter > TheXMLOutputWriter;

  typedef SingletonHolder< XMLFileWriter > TheXMLLogWriter;

} // End namespace Chroma

#endif
