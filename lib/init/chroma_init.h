// -*- C++ -*-
// $Id: chroma_init.h,v 3.0 2006-04-03 04:58:54 edwards Exp $
/*! \file
 *  \brief Initialization of Chroma
 */

#ifndef __CHROMA_INIT_H__
#define __CHROMA_INIT_H__

#include "chromabase.h"

namespace Chroma 
{
  //! Chroma initialisation routine
  void initialize(int* argc, char ***argv);

  //! Chroma finalization routine
  void finalize(void);

  //! Chroma abort routine
  void abort(int i);

  //! Get input file name
  string getXMLInputFileName();

  //! Get output file name
  string getXMLOutputFileName();

  //! Get output logfile name
  string getXMLLogFileName();

  //! Get current working directory
  string getCWD();


  //! Set input file name
  void setXMLInputFileName(const string&);

  //! Set output file name
  void setXMLOutputFileName(const string&);

  //! Set output logfile name
  void setXMLLogFileName(const string&);

  //! Set current working directory
  void setCWD(const string&);


  //! Get xml output instance
  XMLFileWriter& getXMLOutputInstance();
  
  //! Get xml log instance
  XMLFileWriter& getXMLLogInstance();

  /*
  //! Get xml input instance
  XMLReader& getXMLInputInstance();
  */


} // End namespace Chroma

#endif
