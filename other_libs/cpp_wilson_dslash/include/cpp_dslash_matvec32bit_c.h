#ifndef CPP_DSLASH_MATVEC_32BIT_H
#define CPP_DSLASH_MATVEC_32BIT_H

#warning "using C stuff"

#include <cpp_dslash_types.h>

namespace CPlusPlusWilsonDslash {

    using namespace Dslash32BitTypes;

    // Gauge Matrix is 3x3x2 (Natural ordering)
    // HalfSpinor is 3x2x2 (color, spin, reim)
    inline 
    void su3_mult(HalfSpinor res, const GaugeMatrix u, const HalfSpinor src) 
    {

      int re=0;
      int im=1;
 
      // col=0     
      // row = 0
      // spin=0
      res[0][0][re] = u[0][0][re]*src[0][0][re] 
	-u[0][0][im]*src[0][0][im];
      
      res[0][0][im] = u[0][0][re]*src[0][0][im]
	+u[0][0][im]*src[0][0][re];
      
      // spin=1
      res[0][1][re] = u[0][0][re]*src[0][1][re]
	-u[0][0][im]*src[0][1][im];
      
      res[0][1][im] = u[0][0][re]*src[0][1][im]
	+u[0][0][im]*src[0][1][re];
   
      
      // row = 1
      // spin =0
      res[1][0][re] = u[0][1][re]*src[0][0][re]
	-u[0][1][im]*src[0][0][im];
      
      res[1][0][im] = u[0][1][re]*src[0][0][im]
	+u[0][1][im]*src[0][0][re];

      // spin =1
      res[1][1][re] = u[0][1][re]*src[0][1][re]
	-u[0][1][im]*src[0][1][im];
      
      res[1][1][im] = u[0][1][re]*src[0][1][im]
	+u[0][1][im]*src[0][1][re];

      // row = 2
      // spin=0
      res[2][0][re] = u[0][2][re]*src[0][0][re]
	-u[0][2][im]*src[0][0][im];
      
      res[2][0][im] = u[0][2][re]*src[0][0][im]
	+u[0][2][im]*src[0][0][re];
      
      // spin=1
      res[2][1][re] = u[0][2][re]*src[0][1][re]
	-u[0][2][im]*src[0][1][im];
      
      res[2][1][im] = u[0][2][re]*src[0][1][im]
	+u[0][2][im]*src[0][1][re];


      // col=1
      // row=0
      // spin=0
      res[0][0][re] += u[1][0][re]*src[1][0][re]
	-u[1][0][im]*src[1][0][im];
      
      res[0][0][im] += u[1][0][re]*src[1][0][im]
	+u[1][0][im]*src[1][0][re];
      
      // spin=1
      res[0][1][re] += u[1][0][re]*src[1][1][re]
	-u[1][0][im]*src[1][1][im];

      res[0][1][im] += u[1][0][re]*src[1][1][im]
	+u[1][0][im]*src[1][1][re];

      // row=1
      // spin=0
      res[1][0][re] += u[1][1][re]*src[1][0][re]
	-u[1][1][im]*src[1][0][im];
      
      res[1][0][im] += u[1][1][re]*src[1][0][im]
	+u[1][1][im]*src[1][0][re];

      // spin=1
      res[1][1][re] += u[1][1][re]*src[1][1][re]
	-u[1][1][im]*src[1][1][im];
      
      res[1][1][im] += u[1][1][re]*src[1][1][im]
	+u[1][1][im]*src[1][1][re];

      // row=2
      // spin=0
      res[2][0][re] += u[1][2][re]*src[1][0][re]
	-u[1][2][im]*src[1][0][im];
      
      res[2][0][im] += u[1][2][re]*src[1][0][im]
	+u[1][2][im]*src[1][0][re];
      
      // spin=1
      res[2][1][re] += u[1][2][re]*src[1][1][re]
	-u[1][2][im]*src[1][1][im];
      
      res[2][1][im] += u[1][2][re]*src[1][1][im]
	+u[1][2][im]*src[1][1][re];
     


      // col=2
      // row=0
      // spin=0
      res[0][0][re] += u[2][0][re]*src[2][0][re]
	-u[2][0][im]*src[2][0][im];
      
      res[0][0][im] += u[2][0][re]*src[2][0][im]
	+u[2][0][im]*src[2][0][re];
      
      // spin=1
      res[0][1][re] += u[2][0][re]*src[2][1][re]
	-u[2][0][im]*src[2][1][im];
      
      res[0][1][im] += u[2][0][re]*src[2][1][im]
	+u[2][0][im]*src[2][1][re];
      
      // row=1
      // spin=0
      res[1][0][re] += u[2][1][re]*src[2][0][re]
	-u[2][1][im]*src[2][0][im];
      
      res[1][0][im] += u[2][1][re]*src[2][0][im]
	+u[2][1][im]*src[2][0][re];
      
      // spin=1
      res[1][1][re] += u[2][1][re]*src[2][1][re]
	-u[2][1][im]*src[2][1][im];
      
      res[1][1][im] += u[2][1][re]*src[2][1][im]
	+u[2][1][im]*src[2][1][re];
      
      // row = 2
      // spin = 0
      res[2][0][re] += u[2][2][re]*src[2][0][re]
	-u[2][2][im]*src[2][0][im];
      
      res[2][0][im] += u[2][2][re]*src[2][0][im]
	+u[2][2][im]*src[2][0][re];
     
      // spin= 1
      res[2][1][re] += u[2][2][re]*src[2][1][re]
	-u[2][2][im]*src[2][1][im];
      
      res[2][1][im] += u[2][2][re]*src[2][1][im]
	+u[2][2][im]*src[2][1][re];
      
    }

    inline 
      void su3_adj_mult(HalfSpinor res, GaugeMatrix u, HalfSpinor src)
    {
      int re=0;
      int im=1;

      // col = 0
      // row =0
      // spin =0 
      res[0][0][re] = u[0][0][re]*src[0][0][re]
	+ u[0][0][im]*src[0][0][im];
      res[0][0][im] = u[0][0][re]*src[0][0][im]
	- u[0][0][im]*src[0][0][re];
      
      // spin=1
      res[0][1][re] = u[0][0][re]*src[0][1][re]
	+ u[0][0][im]*src[0][1][im];
      res[0][1][im] = u[0][0][re]*src[0][1][im]
	  - u[0][0][im]*src[0][1][re];
      
      // col = 1
      // spin =0 
      res[0][0][re] += u[0][1][re]*src[1][0][re]
	+ u[0][1][im]*src[1][0][im];
      
      res[0][0][im] += u[0][1][re]*src[1][0][im]
	- u[0][1][im]*src[1][0][re];
      
      // spin =1
      res[0][1][re] += u[0][1][re]*src[1][1][re]
	+ u[0][1][im]*src[1][1][im];
      
      res[0][1][im] += u[0][1][re]*src[1][1][im]
	- u[0][1][im]*src[1][1][re];

      // col=2
      // spin =0 
      res[0][0][re] += u[0][2][re]*src[2][0][re]
	+ u[0][2][im]*src[2][0][im];
      
      res[0][0][im] += u[0][2][re]*src[2][0][im]
	- u[0][2][im]*src[2][0][re];
      
      // spin =1
      res[0][1][re] += u[0][2][re]*src[2][1][re]
	+ u[0][2][im]*src[2][1][im];
      
      res[0][1][im] += u[0][2][re]*src[2][1][im]
	- u[0][2][im]*src[2][1][re];



      // row =1 
      // col =0
      // spin =0 
      res[1][0][re] = u[1][0][re]*src[0][0][re]
	+ u[1][0][im]*src[0][0][im];
      res[1][0][im] = u[1][0][re]*src[0][0][im]
	- u[1][0][im]*src[0][0][re];
	
      // spin=1
      res[1][1][re] = u[1][0][re]*src[0][1][re]
	+ u[1][0][im]*src[0][1][im];
      res[1][1][im] = u[1][0][re]*src[0][1][im]
	- u[1][0][im]*src[0][1][re];
      
      // col=1
      // spin =0 
      res[1][0][re] += u[1][1][re]*src[1][0][re]
	+ u[1][1][im]*src[1][0][im];
      
      res[1][0][im] += u[1][1][re]*src[1][0][im]
	- u[1][1][im]*src[1][0][re];
      
      // spin =1
      res[1][1][re] += u[1][1][re]*src[1][1][re]
	+ u[1][1][im]*src[1][1][im];
      
      res[1][1][im] += u[1][1][re]*src[1][1][im]
	- u[1][1][im]*src[1][1][re];

      // col=2 
      // spin =0 
      res[1][0][re] += u[1][2][re]*src[2][0][re]
	+ u[1][2][im]*src[2][0][im];
      
      res[1][0][im] += u[1][2][re]*src[2][0][im]
	- u[1][2][im]*src[2][0][re];
      
      // spin =1
      res[1][1][re] += u[1][2][re]*src[2][1][re]
	+ u[1][2][im]*src[2][1][im];
      
      res[1][1][im] += u[1][2][re]*src[2][1][im]
	- u[1][2][im]*src[2][1][re];


      // row=2
      // spin =0 
      res[2][0][re] = u[2][0][re]*src[0][0][re]
	+ u[2][0][im]*src[0][0][im];
      res[2][0][im] = u[2][0][re]*src[0][0][im]
	- u[2][0][im]*src[0][0][re];
	
      // spin=1
      res[2][1][re] = u[2][0][re]*src[0][1][re]
	+ u[2][0][im]*src[0][1][im];
      res[2][1][im] = u[2][0][re]*src[0][1][im]
	- u[2][0][im]*src[0][1][re];
      
      // col=1
      // spin =0 
      res[2][0][re] += u[2][1][re]*src[1][0][re]
	+ u[2][1][im]*src[1][0][im];
	
      res[2][0][im] += u[2][1][re]*src[1][0][im]
	- u[2][1][im]*src[1][0][re];
      
      // spin =1
      res[2][1][re] += u[2][1][re]*src[1][1][re]
	+ u[2][1][im]*src[1][1][im];
      
      res[2][1][im] += u[2][1][re]*src[1][1][im]
	- u[2][1][im]*src[1][1][re];
      
      // col=2
      // spin =0 
      res[2][0][re] += u[2][2][re]*src[2][0][re]
	+ u[2][2][im]*src[2][0][im];
      
      res[2][0][im] += u[2][2][re]*src[2][0][im]
	- u[2][2][im]*src[2][0][re];
	
      // spin =1
      res[2][1][re] += u[2][2][re]*src[2][1][re]
	+ u[2][2][im]*src[2][1][im];
      
      res[2][1][im] += u[2][2][re]*src[2][1][im]
	- u[2][2][im]*src[2][1][re];

    }
  
}

#endif

    



