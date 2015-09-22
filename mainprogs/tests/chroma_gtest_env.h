#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstdio>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "chroma.h"
#include "gtest/gtest.h"

using namespace Chroma;

/* Set up Chroma */
class ChromaEnvironment : public ::testing::Environment {
public: 
  ChromaEnvironment(int* argc, char ***argv) {
    Chroma::initialize(argc, argv);
    QDPIO::cout << "Linkage = " << linkageHack() << std::endl;
  }
  
  ~ChromaEnvironment() {
    Chroma::finalize();
 }

  bool linkageHack(void)
  {
    bool foo = true;
    
    // Inline Measurements
    foo &= InlineAggregateEnv::registerAll();
    foo &= GaugeInitEnv::registerAll();
    
    return foo;
  }
  
};

