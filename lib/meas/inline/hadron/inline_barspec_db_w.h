// -*- C++ -*-
/*! \file
 * \brief Inline hadron spectrum calculations
 *
 * Hadron spectrum calculations
 */

#ifndef __inline_barspec_h__
#define __inline_barspec_h__

#include <map>
#include <vector>

#include "chromabase.h"
//#include "../utils/gammaRotations.h"
#include "util/ferm/diractodr.h"

#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineBarSpecEnv 
  {
    bool registerAll();


    //! Parameter structure
    /*! \ingroup inlinehadron */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long frequency;

      struct SpinTerms_t{
	multi1d<int> spin ;
	double weight ;
      };

      struct SpinWF_t{
	double norm ;
	multi1d<SpinTerms_t> terms ;
      } ;


      struct Operators_t{
	string name;

	SpinWF_t spinWF ;
      };

      struct State_t{
	string name;
	multi1d<string> flavor ;
	multi1d<Operators_t> ops ;        // holds the operators
	int spin ;
	string db ;
      };
      struct Param_t
      {
	//bool MesonP;             // Meson spectroscopy
	//bool CurrentP;           // Meson currents
	//bool BaryonP;            // Baryons spectroscopy

	bool time_rev;           // Use time reversal in baryon spectroscopy

	int mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.
	bool avg_equiv_mom;      // average over equivalent momenta

	multi1d<State_t> states ;        // holds the states
      
	string ensemble ; // a string describing this ensemble 
      } param;

      struct NamedObject_t
      {
	std::string  gauge_id;           /*!< Input gauge field */

	struct Props_t
	{
	  std::string  up_id;
	  std::string  down_id;
	  std::string  strange_id;
	  std::string  charm_id;
	};

	multi1d<Props_t> props;

      } named_obj;

      std::string xml_file;  // Alternate XML file pattern

    };


    namespace BarSpec
    {
      vector<int> permutation(int k, const vector<int>& s) ;
      int permutation_sign(int k, int n) ;

      class SpinWF_t{
      private:
      
	int factorial(int n){
	  int f = 1 ;
	  for(int i(1);i<=n;i++)
	    f *= i ;
	  return f ;
	}


      public:
	double norm ;
	vector<Params::SpinTerms_t> terms ;
      
	SpinWF_t(const Params::SpinWF_t& ss):norm(ss.norm){
	  for(int t(0);t<ss.terms.size();t++)
	    terms.push_back(ss.terms[t]) ;
	  //make spin components zero based
	  for(int t(0);t<terms.size();t++){
	    for(int s(0);s<terms[t].spin.size();s++)
	      terms[t].spin[s]--;
	  }
	}

	void permutations(const multi1d<string>& flavor){

	  //map<string,int> f_count;
	  map<string,vector<int> > f_pos; // positions of active flavors
	  //find the common flavors
	  // f_pos  : list of positions on which each flavor appears 
	  // labels : holds the list of unique flavors in the operator
	  vector<string> labels; 
	  for(int f(0);f<flavor.size();f++){
	    if(f_pos.find(flavor[f]) == f_pos.end()){
	      labels.push_back(flavor[f]);
	    }
	    f_pos[flavor[f]].push_back(f);
	  }

	  vector<Params::SpinTerms_t> perm_terms ;
	  for(int f(0);f<labels.size();f++){
	    int Nterms = terms.size();
	    //vector<SpinTerms_t> tmp_terms = terms ;
	    // perform up qark permutations of each term	  
	    int Nperms = factorial(f_pos[labels[f]].size()) ;
	    QDPIO::cout<<"   SpinWF_t::"<<__func__
		       <<": Flavor "<<labels[f]<<" needs "<<Nperms
		       <<": permutations"<<endl ;
	    if(Nperms>1){
	      for(int t(0);t<Nterms;t++){
		Params::SpinTerms_t foo = terms[t] ;
		vector<int> spin ;
		for(int i(0);i<f_pos[labels[f]].size();i++)
		  spin.push_back(terms[t].spin[f_pos[labels[f]][i]]);
		//QDPIO::cout<<"Original : " ;
		//for(int m(0);m<spin.size();m++) 
		//QDPIO::cout<<spin[m]<<" " ;
		//QDPIO::cout<<endl ;
		for(int k(0);k<Nperms;k++){
		  vector<int> perm = permutation(k,spin);
		  //for(int m(0);m<perm.size();m++)
		  //  QDPIO::cout<<perm[m]<<" " ;
		  //QDPIO::cout<<endl ;

		  for(int i(0);i<f_pos[labels[f]].size();i++)
		    foo.spin[f_pos[labels[f]][i]] = perm[i] ;
		  perm_terms.push_back(foo) ;
		}//loop over permutations
	      }// loop over terms 
	    }// If Nperms>1
	  }// loop over flavors 
	  terms = perm_terms ;
	  // Check what we did ...
	  for(int t(0);t<terms.size();t++){
	    QDPIO::cout<<"   SpinWF_t::"<<__func__
		       <<": weight: "<<terms[t].weight<<" | " ;
	    QDPIO::cout<<"spin  : " ;
	    for(int m(0);m<terms[t].spin.size();m++)
	      QDPIO::cout<<flavor[m]<<"("<<terms[t].spin[m]<<") " ;
	    QDPIO::cout<<">"<<endl ;

	  }
	  //QDPIO::cout<<"Done with permutations"<<endl ;
	}
      
      } ;
    
      //My own propagator class                                              
      template<int N>
      class BasePropagator
      {
      public:
	multi2d<multi2d<LatticeComplex> > p ;
      
	BasePropagator(){
	  p.resize(N,N);
	  for(int s(0);s<N;s++)
	    for(int ss(0);ss<N;ss++)
	      p(s,ss).resize(Nc,Nc) ;
	}
      
	void ConvertProp(const LatticePropagator& prop, int spin_offset){
	  QDPIO::cout<<__func__<<": Converting to Dirac Basis"<<endl ;
	  SpinMatrix U = DiracToDRMat();
	  LatticePropagator foo = adj(U)*prop*U ;
	
	  p.resize(N,N);
	  int Noff = N+spin_offset;
	  if(Noff>Ns)
	    Noff=Ns ; 
	  for(int s(spin_offset);s<Noff;s++)
	    for(int ss(0);ss<Noff;ss++){
	      p(s,ss).resize(Nc,Nc) ;
	      LatticeColorMatrix Cmat = peekSpin(foo,s,ss) ;
	      for(int c(0);c<Nc;c++)
		for(int cc(0);cc<Nc;cc++)
		  p(s,ss)(c,cc) = peekColor(Cmat,c,cc);
	    }
	
	}
      
	void ConvertProp(const LatticePropagator& prop){
	  ConvertProp(prop,0);
	}

	void ConvertPropUpper(const LatticePropagator& prop){
	  ConvertProp(prop,0);
	}

	void ConvertPropLowerr(const LatticePropagator& prop){
	  ConvertProp(prop,2);
	}

      };

      typedef BasePropagator<Ns  >  RPropagator ;
      typedef BasePropagator<Ns/2> NRPropagator ;


      void contract(LatticeComplex& latC, 
		    const RPropagator& q1,
		    const RPropagator& q2,
		    const RPropagator& q3,
		    const SpinWF_t& snk,
		    const SpinWF_t& src);

		    


    }// namespace BarSpec


    //! Inline measurement of hadron spectrum
    /*! \ingroup inlinehadron */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the measurement
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    protected:
      //! Do the measurement
      void func(const unsigned long update_no,
		XMLWriter& xml_out); 

    private:
      Params params;
    };

  } // namespace InlineBarSpecEnv

} // namespace Chroma

#endif
