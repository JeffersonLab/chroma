// $Id: time_dslash_3d.cc,v 1.1 2008-08-26 13:24:29 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "unittest.h"
#include "testvol.h"

#include "timeDslash3D.h"

using namespace QDP;

int main(int argc, char **argv)
{
  // Initialize QDP++ with argc, and argv. Set Lattice Dimensions

  // Initialize UnitTest jig
  TestRunner  tests(&argc, &argv, nrow_in);
  tests.addTest(new timeDslash3D(), "timeDslash3D" );

  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();

}

