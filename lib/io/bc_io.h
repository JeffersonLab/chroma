#ifndef BC_IO_H
#define BC_IO_H
#include "chromabase.h"

using namespace QDP;

// Commonly used Boundaries for fermions.
enum SimpleBCType_t { BC_ANTIPERIODIC = -1,  // Antiperiodic = -1
		      BC_DIRICHLET,          // Dirichlet    =  0
		      BC_PERIODIC };         // Periodic     =  1

void read(XMLReader& xml, const std::string& path, SimpleBCType_t& t);
void write(XMLWriter& xml, const std::string& path, const SimpleBCType_t& t);


#endif
