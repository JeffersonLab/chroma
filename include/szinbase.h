// $Id: szinbase.h,v 1.1 2003-02-15 05:57:00 edwards Exp $
//
// Absolute basic stuff to use szin

#include "qdp.h"
#include "linearop.h"

using namespace QDP;

//void CREATE_NAMELIST(const char *s, ...);
//inline void WRITE_NAMELIST(...) {}
#define CREATE_NAMELIST(a,b) 

// #define WRITE_NAMELIST(grp,...) Write(cout, #grp , ##__VA_ARGS__)


const float fuzz = 1.0e-5;
const float twopi = 6.283185307179586476925286;

#define START_CODE(a)
#define END_CODE(a)

#define TO_REAL(a) float(a)
