// -*- C++ -*-
// $Id: xmllog_io.h,v 1.4 2005-03-02 00:42:37 edwards Exp $
/*! \file
 * \brief Singleton instances of xml output
 */

#ifndef _XMLLOG_IO_H_
#define _XMLLOG_IO_H_

#include "chromabase.h"
#include "singleton.h"


namespace Chroma 
{

  typedef SingletonHolder< XMLFileWriter > TheXMLOutputWriter;

  typedef SingletonHolder< XMLFileWriter > TheXMLLogWriter;

} // End namespace Chroma

#endif
