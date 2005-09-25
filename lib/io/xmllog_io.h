// -*- C++ -*-
// $Id: xmllog_io.h,v 2.0 2005-09-25 21:04:32 edwards Exp $
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
