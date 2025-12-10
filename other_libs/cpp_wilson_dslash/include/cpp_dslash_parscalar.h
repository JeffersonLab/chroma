#ifndef CPP_DSLASH_PARSCALAR_H
#define CPP_DSLASH_PARSCALAR_H

#include <cstddef>


/* This is the scalar dslash */ 
/* Use the scalar types */
#include <cpp_dslash_types.h>
#include <shift_table_parscalar.h>
#include <shift_table_3d_parscalar.h>
#include <tables_parscalar.h>
#include <qmp.h>

namespace CPlusPlusWilsonDslash {

  /* Forward Declare the shift table */
  // class ShiftTable;


  /* Specialize - Dslash of float */
  template<typename T> 
    class Dslash {};

  template<typename T> 
    class Dslash3D {};


  template<>
    class Dslash<float> {
    public:
      
      /* Use the 32 bit types here */
      typedef Dslash32BitTypes::GaugeMatrix GaugeMatrix;   /* 3x3 complex matrices */
      typedef Dslash32BitTypes::FourSpinor FourSpinor;    /* 4x3 dirac spinors */
      typedef Dslash32BitTypes::HalfSpinor HalfSpinor;
      
      /* Constructor */
      /* Parameters: latt_size[] - total lattice size */
      /*             getSiteCoords() - geometry function from external package */
      /*             getLinearSiteIndex() - geometry function from external package */
      /*             nodeNum() - geometry function from external package */
      
      Dslash(const int latt_size[],      
	     void (*getSiteCoords)(int coord[], int node, int linearsite),
	     int (*getLinearSiteIndex)(const int coord[]),
	     int (*nodeNum)(const int coord[])
	     );
      
      /* Destructor */
      ~Dslash();
      
      /* Apply the operator */
      void operator()(float* res, 
		      float* psi, 
		      float* u, 
		      int isign,
		      int cb);
      
    private:
      // ShiftTable* s;             // Some Offset Table
      
      // Temporary space...
      static ShiftTable<HalfSpinor>* s_tab;
      static DslashTables<HalfSpinor,4>* tab;
      
      // Tab
      Dslash();
    }; // Class
  

/* Specialize - Dslash of double */
  template<>
    class Dslash<double> {
    public:
      
      typedef Dslash64BitTypes::GaugeMatrix GaugeMatrix;   /* 3x3 complex matrices */
      typedef Dslash64BitTypes::FourSpinor FourSpinor;    /* 4x3 dirac spinors */
      typedef Dslash64BitTypes::HalfSpinor HalfSpinor;
      
      /* Constructor */
      /* Parameters: latt_size[] - total lattice size */
      /*             getSiteCoords() - geometry function from external package */
      /*             getLinearSiteIndex() - geometry function from external package */
      /*             nodeNum() - geometry function from external package */
      
      Dslash(const int latt_size[],      
	     void (*getSiteCoords)(int coord[], int node, int linearsite),
	     int (*getLinearSiteIndex)(const int coord[]),
	     int (*nodeNum)(const int coord[])
	     );
      
      /* Destructor */
      ~Dslash();
      
      /* Apply the operator */
      void operator()(double* res, 
		      double* psi, 
		      double* u, 
		      int isign,
		      int cb);
      
    private:
      static ShiftTable<HalfSpinor>* s_tab;
      static DslashTables<HalfSpinor,4>* tab;
      
      // Tab
      Dslash();
    }; // Class






  template<>
    class Dslash3D<float> {
    public:
      
      /* Use the 32 bit types here */
      typedef Dslash32BitTypes::GaugeMatrix GaugeMatrix;   /* 3x3 complex matrices */
      typedef Dslash32BitTypes::FourSpinor FourSpinor;    /* 4x3 dirac spinors */
      typedef Dslash32BitTypes::HalfSpinor HalfSpinor;
      
      /* Constructor */
      /* Parameters: latt_size[] - total lattice size */
      /*             getSiteCoords() - geometry function from external package */
      /*             getLinearSiteIndex() - geometry function from external package */
      /*             nodeNum() - geometry function from external package */
      
      Dslash3D(const int latt_size[],      
	     void (*getSiteCoords)(int coord[], int node, int linearsite),
	     int (*getLinearSiteIndex)(const int coord[]),
	     int (*nodeNum)(const int coord[])
	     );
      
      /* Destructor */
      ~Dslash3D();
      
      /* Apply the operator */
      void operator()(float* res, 
		      float* psi, 
		      float* u, 
		      int isign,
		      int cb);
      
    private:
      // Temporary space...
     static ShiftTable3D<HalfSpinor>* s_tab;
     static DslashTables<HalfSpinor,3>* tab;
      
      // Tab
      Dslash3D();
    }; // Class
  

/* Specialize - Dslash of double */
  template<>
    class Dslash3D<double> {
    public:
      
      typedef Dslash64BitTypes::GaugeMatrix GaugeMatrix;   /* 3x3 complex matrices */
      typedef Dslash64BitTypes::FourSpinor FourSpinor;    /* 4x3 dirac spinors */
      typedef Dslash64BitTypes::HalfSpinor HalfSpinor;
      
      /* Constructor */
      /* Parameters: latt_size[] - total lattice size */
      /*             getSiteCoords() - geometry function from external package */
      /*             getLinearSiteIndex() - geometry function from external package */
      /*             nodeNum() - geometry function from external package */
      
      Dslash3D(const int latt_size[],      
	     void (*getSiteCoords)(int coord[], int node, int linearsite),
	     int (*getLinearSiteIndex)(const int coord[]),
	     int (*nodeNum)(const int coord[])
	     );
      
      /* Destructor */
      ~Dslash3D();
      
      /* Apply the operator */
      void operator()(double* res, 
		      double* psi, 
		      double* u, 
		      int isign,
		      int cb);
      
    private:
      static ShiftTable3D<HalfSpinor>* s_tab;
      static DslashTables<HalfSpinor,3>* tab;
      
      // Tab
      Dslash3D();
    }; // Class


} // Namespace

#endif
