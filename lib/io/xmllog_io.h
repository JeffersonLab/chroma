// -*- C++ -*-
// $Id: xmllog_io.h,v 1.5 2005-03-02 18:32:05 bjoo Exp $
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
