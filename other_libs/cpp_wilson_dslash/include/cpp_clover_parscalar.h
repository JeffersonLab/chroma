#ifndef CPP_CLOVER_PARSCALAR_H
#define CPP_CLOVER_PARSCALAR_H

#include <cstddef>


/* This is the scalar dslash */ 
/* Use the scalar types */
#include <cpp_dslash_types.h>
#include <cpp_clover_types.h>

#include <shift_table_parscalar.h>
#include <tables_parscalar.h>
#include <qmp.h>

using namespace CPlusPlusWilsonDslash;
namespace CPlusPlusClover {

  // Never Created
  template<typename T> 
    class CloverSchur4D {};
  
  // Single precision
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
      void operator()(float* res, 
		      const float* psi, 
		      const float* u, 
		      const float* clov_oo,
		      const float* invclov_ee,
		      int isign);
      
    private:
      // Temporary space...
      ShiftTable<HalfSpinor>* s_tab;
      DslashTables<HalfSpinor,4>* tab;
      FourSpinor* xt_spinor; // Unaligned
      FourSpinor* t_spinor; //  Aligned (Use this to store A^{-1}_ee D_oe psi)
      FourSpinor* xt_spinor2; // Unaligned
      FourSpinor* t_spinor2;  // Aligned (Use this to store A_oo psi)

      // Tab
      CloverSchur4D();
    }; // Class
  

  
  // Single precision
  template<>
    class CloverSchur4D<double> {
    public:
      
      /* Use the 64 bit types here */
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
      void operator()(double* res, 
		      const double* psi, 
		      const double* u, 
		      const double* clov_oo,
		      const double* invclov_ee,
		      int isign);
      
    private:
      // Temporary space...
      ShiftTable<HalfSpinor>* s_tab;
      DslashTables<HalfSpinor,4>* tab;
      FourSpinor* xt_spinor; // Unaligned
      FourSpinor* t_spinor; //  Aligned (Use this to store A^{-1}_ee D_oe psi)
      FourSpinor* xt_spinor2; // Unaligned
      FourSpinor* t_spinor2;  // Aligned (Use this to store A_oo psi)

      // Tab
      CloverSchur4D();
    }; // Class
  


} // Namespace

#endif
