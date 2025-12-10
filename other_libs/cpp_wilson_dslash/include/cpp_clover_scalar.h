#ifndef CPP_CLOVER_SCALAR_H
#define CPP_CLOVER_SCALAR_H

#include <cstddef>


/* This is the scalar dslash */ 
/* Use the scalar types */
#include <cpp_dslash_types.h>
#include <cpp_clover_types.h>
#include <cpp_dslash_scalar.h>

using namespace CPlusPlusWilsonDslash;

namespace CPlusPlusClover {

  /* Forward Declare the shift table */
  
  
  /* This we will never use */
  template<typename FT>
  class CloverSchur4D { };  /* Empty class */

  /* Specialize - a Single Prec clover op */
  template<>  
  class CloverSchur4D<float> {
  public:

    /* Use the 32 bit types here */
    typedef Dslash32BitTypes::GaugeMatrix GaugeMatrix;   /* 3x3 complex matrices */
    typedef Dslash32BitTypes::FourSpinor FourSpinor;    /* 4x3 dirac spinors */
    typedef Dslash32BitTypes::HalfSpinor HalfSpinor; 
    typedef Clover32BitTypes::CloverTerm CloverTerm;

    /* Constructor */
    /* Parameters: latt_size[] - total lattice size */
    /*             getSiteCoords() - geometry function from external package */
    /*             getLinearSiteIndex() - geometry function from external package */
    /*             nodeNum() - geometry function from external package */

    CloverSchur4D(const int latt_size[],      
		  void (*getSiteCoords)(int coord[], int node, int linearsite),
		  int (*getLinearSiteIndex)(const int coord[]),
		  int (*nodeNum)(const int coord[])
		  );
  
    /* Destructor */
    ~CloverSchur4D();
  
    /* Apply the operator */
    void operator()(float* res,  // Result 
		    const float* psi,  // Src
		    const float* u,    // Gauge Field
		    const float* clov_oo, // Clover Term on Odd Sites
		    const float* invclov_ee, // Clover Inverse on Even Sites
		    int isign);

  private:
    CPlusPlusWilsonDslash::ShiftTable* s;
    FourSpinor *xt_spinor;
    FourSpinor *t_spinor;

    // Hide Free Constructor 
    CloverSchur4D();
  }; // Class

  /* Specialize - a Double Prec clover op */
  template<>  
  class CloverSchur4D<double> {
  public:

    /* Use the 32 bit types here */
    typedef Dslash64BitTypes::GaugeMatrix GaugeMatrix;   /* 3x3 complex matrices */
    typedef Dslash64BitTypes::FourSpinor FourSpinor;    /* 4x3 dirac spinors */
    typedef Dslash64BitTypes::HalfSpinor HalfSpinor; 
    typedef Clover64BitTypes::CloverTerm CloverTerm;

    /* Constructor */
    /* Parameters: latt_size[] - total lattice size */
    /*             getSiteCoords() - geometry function from external package */
    /*             getLinearSiteIndex() - geometry function from external package */
    /*             nodeNum() - geometry function from external package */

    CloverSchur4D(const int latt_size[],      
		  void (*getSiteCoords)(int coord[], int node, int linearsite),
		  int (*getLinearSiteIndex)(const int coord[]),
		  int (*nodeNum)(const int coord[])
		  );
  
    /* Destructor */
    ~CloverSchur4D();
  
    /* Apply the operator */
    void operator()(double* res,  // Result 
		    const double* psi,  // Src
		    const double* u,    // Gauge Field
		    const double* clov_oo, // Clover Term on Odd Sites
		    const double* invclov_ee, // Clover Inverse on Even Sites
		    int isign);

  private:
    CPlusPlusWilsonDslash::ShiftTable* s;
    FourSpinor *xt_spinor;
    FourSpinor *t_spinor;

    // Hide Free Constructor 
    CloverSchur4D();

  }; // Class


} // Namespace

#endif
