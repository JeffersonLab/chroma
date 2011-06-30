#include "chromabase.h"
#include "update/molecdyn/integrator/integrator_shared.h"

using namespace QDP;
namespace Chroma {
  namespace IntegratorShared { 

    AbsComponentIntegrator< multi1d<LatticeColorMatrix>,
			    multi1d<LatticeColorMatrix> >* 
    createSubIntegrator(const std::string& subintegrator_xml) 
    {
      
      std::istringstream is( subintegrator_xml );
      XMLReader top(is);
      
      std::string subint_name;
      try { 
	read(top, "/SubIntegrator/Name", subint_name);
      }
      catch( const std::string& e ) {
	QDPIO::cerr << "Failed to extract name of subintegrator in LatColMatSTSLeapfrogRecursiveIntegrator: " << e << endl;
	QDP_abort(1);
      }
      std::string root="/SubIntegrator";
    
      AbsComponentIntegrator< multi1d<LatticeColorMatrix>,
	multi1d<LatticeColorMatrix> >* ret_val=
	TheMDComponentIntegratorFactory::Instance().createObject(subint_name, top, root);
      return ret_val;
    }

    void bindMonomials(const multi1d<std::string>& monomial_ids,
		       multi1d< MonomialPair >& monomials ) { 
      QDPIO::cout << "Binding Monomials" << endl;
      QDPIO::cout << "There are " << monomial_ids.size() << " IDs to bind " << endl; 

      monomials.resize(monomial_ids.size());
      for(int i=0; i < monomial_ids.size(); i++) { 
	try { 

	  monomials[i].mon=
	    TheNamedObjMap::Instance().getData< Handle< Monomial< multi1d<LatticeColorMatrix> , multi1d<LatticeColorMatrix> > > >(monomial_ids[i]);

	  monomials[i].id = monomial_ids[i];

	  QDPIO::cout << "Monomial with ID = " << monomial_ids[i] << " bound" << endl;
	}
	catch(const std::string& e) { 
	  QDPIO::cout << "Caught exception with message : " << e << endl;
	  QDP_abort(1);
	}
	catch(std::bad_cast) { 
	  QDPIO::cout << "Failed to cast to return type in bind_monomials" << endl;
	  QDP_abort(1);
	}

      } // For loop 
      QDPIO::cout << "All monomials successfully bound" << endl;
    }


  }

}
