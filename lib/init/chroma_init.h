// $Id: chroma_init.h,v 1.3 2005-01-17 21:08:29 edwards Exp $

#ifndef __CHROMA_INIT_H__
#define __CHROMA_INIT_H__

#include "chromabase.h"

namespace Chroma {

  //! Chroma initialisation routine
  void ChromaInitialize(int* argc, char ***argv);

  //! Chroma finalization routine
  void ChromaFinalize(void);

  //! Chroma abort routine
  void ChromaAbort(int i);

}; // End namespace Chroma

#endif
