/*! \file
 *  \brief NEF fermion action parameters
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/nef_fermact_params_w.h"

#include "io/param_io.h"

namespace Chroma
{

  //! Default constructor
  NEFFermActParams::NEFFermActParams() : b5(1), c5(1)
  {	  
    Mass = Real(0);
    OverMass= Real(0);
    N5=1;
    b5[0]=Real(1);
    c5[0]=Real(1);
  }

  //! Read parameters
  NEFFermActParams::NEFFermActParams(XMLReader& xml, const std::string& path)
  {	  
    XMLReader paramtop(xml, path);
    read(paramtop, "Mass", Mass);
    read(paramtop, "OverMass", OverMass);
    read(paramtop, "N5", N5);
    read(paramtop, "b5", b5);
    read(paramtop, "c5", c5);

    // For the regular Moebius there is 1 value of b5 and one of c5
    // For the general one there is 1 per 5th-dim slice.
    // Here I check that either the array is length N5
    // or if it is length 1 (regular Moebius) I replicate 
    // the one value N5 times.

    if ( b5.size() != N5 ) { 
     
      if ( b5.size() == 1 ) { // Only one value entered
	Real b5_f = b5[0]; // get the value
	b5.resize(N5);     // resize array 
	for(int s=0; s < N5; s++) {  // replicate
	  b5[s] = b5_f;
	}
      }
      else { 
	// Not length=N5 and not length = 1: Bomb
	QDPIO::cerr << "b5 must have either lenght 1 or N5" << std::endl;
	QDP_abort(1);
      }

    }

    // Repeat for c5
    if ( c5.size() != N5 ) { 
     
      if ( c5.size() == 1 ) { // Only one value entered
	Real c5_f = c5[0]; // get the value
	c5.resize(N5);     // resize array 
	for(int s=0; s < N5; s++) {  // replicate
	  c5[s] = c5_f;
	}
      }
      else { 
	// Not length=N5 and not length = 1: Bomb
	QDPIO::cerr << "b5 must have either lenght 1 or N5" << std::endl;
	QDP_abort(1);
      }
    }

  }

  //! Read parameters
  void read(XMLReader& xml, const std::string& path, NEFFermActParams& param)
  {	  
    NEFFermActParams tmp(xml, path);	
    param = tmp;
  }


  //! Write parameters
  void write(XMLWriter& xml, const std::string& path, const NEFFermActParams& param)
  {
    push(xml, path);

    write(xml, "Mass", param.Mass);
    write(xml, "OverMass", param.OverMass);
    write(xml, "N5", param.N5);
    write(xml, "b5", param.b5);
    write(xml, "c5", param.c5);

    pop(xml);
  }


}

