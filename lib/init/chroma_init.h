// -*- C++ -*-
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
  std::string getXMLInputFileName();

  //! Get output file name
  std::string getXMLOutputFileName();

  //! Get output logfile name
  std::string getXMLLogFileName();

  //! Get current working directory
  std::string getCWD();


  //! Set input file name
  void setXMLInputFileName(const std::string&);

  //! Set output file name
  void setXMLOutputFileName(const std::string&);

  //! Set output logfile name
  void setXMLLogFileName(const std::string&);

  //! Set current working directory
  void setCWD(const std::string&);


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
