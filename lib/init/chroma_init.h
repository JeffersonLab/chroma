#ifndef __CHROMA_INIT_H__
#define __CHROMA_INIT_H__

#include "chroma.h"
#include <string>


namespace Chroma {

  //! Chroma initialisation routine
  void ChromaInitialize(int* argc, char ***argv);

  //! Chroma finalization routine
  void ChromaFinalize(void);

  //! Chroma abort routine
  void ChromaAbort(int i);

}; // End namespace Chroma

#endif
