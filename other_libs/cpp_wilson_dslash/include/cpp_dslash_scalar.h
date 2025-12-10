#ifndef CPP_DSLASH_SCALAR_H
#define CPP_DSLASH_SCALAR_H

#include <cstddef>


/* This is the scalar dslash */ 
/* Use the scalar types */
#include <cpp_dslash_types.h>

namespace CPlusPlusWilsonDslash {

  /* Forward Declare the shift table */
  class ShiftTable;
  
  /* This we will never use */
  template<typename FT>
  class Dslash { };  /* Empty class */

  template<typename FT>
    class Dslash3D {};

  /* Specialize - Dslash of float */
  template<>  /* The teplate is the floating point type */
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

    //   int getPathSite(int site) const;

  private:
    ShiftTable* s;

    // Hide Free Constructor 
    Dslash();
  }; // Class


  template<>  /* The teplate is the floating point type */
  class Dslash<double> {
  public:

    /* Use the 32 bit types here */
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


    //   int getPathSite(int site) const;

    
  private:
    ShiftTable* s;

    // Hide Free Constructor 
    Dslash();
  }; // Class


  /* Forward declare this */
  class ShiftTable3D;


  /* Specialize - Dslash3D of float */
  template<>  /* The teplate is the floating point type */
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
    ShiftTable3D* s;
    
    // Hide Free Constructor 
    Dslash3D();
  }; // Class


  template<>  /* The teplate is the floating point type */
  class Dslash3D<double> {
  public:

    /* Use the 32 bit types here */
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
    ShiftTable3D* s;

    // Hide Free Constructor 
    Dslash3D();
  }; // Class
} // Namespace

#endif
