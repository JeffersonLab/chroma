// -*- C++ -*-
// $Id: exact_hamiltonian.h,v 3.5 2008-05-21 17:07:50 bjoo Exp $
/*! \file
 * \brief Exact Hamiltonians
 */

#ifndef EXACT_HAMILTONIAN_H
#define EXACT_HAMILTONIAN_H

#include "chromabase.h"
#include "handle.h"

#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"

#include "io/xmllog_io.h"
#include "io/monomial_io.h"
#include "meas/inline/io/named_objmap.h"

//PSILVA
#include "tower_array.h"
#include "pq_traits.h"

namespace Chroma 
{



  //! Parameter structure for new Hamiltonian
  /*! @ingroup hamilton */
  struct ExactHamiltonianParams { 
    
    //! Constructor 
    ExactHamiltonianParams(XMLReader& xml, const std::string& path); 
    multi1d<std::string> monomial_ids; /*!< list of monomial IDs */

  };

  //! Read the parameters for the Hamiltonian
  void read(XMLReader& xml, const std::string& path, ExactHamiltonianParams& p);

  //! Write the parameters for the Hamiltonian
  void write(XMLWriter& xml, const std::string& path, const ExactHamiltonianParams& p);


  //! The Exact Hamiltonian Class - supplies internal field refreshment and energy calculations
  /*! @ingroup hamilton */
  class ExactHamiltonian : public AbsHamiltonian< multi1d<LatticeColorMatrix>, 
		      multi1d<LatticeColorMatrix> >
  {
  public:

    //! Construct from a list of string monomial_ids
    ExactHamiltonian(const multi1d<std::string>& monomial_ids_)  {
      create(monomial_ids_);
    }
   
    //! Construct from a parameter structure
    ExactHamiltonian(const ExactHamiltonianParams& p) {
      create(p.monomial_ids);
    }

    //! Copy constructor
    ExactHamiltonian(const ExactHamiltonian& H) : monomials(H.monomials) {}

    //! Destructor 
    ~ExactHamiltonian(void) {}

    //! Internal Field Refreshment 
    void refreshInternalFields(const AbsFieldState<multi1d<LatticeColorMatrix>,multi1d<LatticeColorMatrix> >& s)  
    { 
      START_CODE();
      for(int i=0; i < monomials.size(); i++) {
	monomials[i]->refreshInternalFields(s);
      }
      END_CODE();
    }
 

    Double mesKE(const AbsFieldState< 
	       multi1d<LatticeColorMatrix>, 
	       multi1d<LatticeColorMatrix> > &s
	       ) const
    {
      START_CODE();

      /* Accumulate KE per site */
      multi1d<LatticeDouble> ke_per_site(Nd);

      
      /* Now add on the local Norm2 of the momenta for each link */
      for(int mu=0; mu < Nd; mu++) { 
	ke_per_site[mu] = -Double(4);
	ke_per_site[mu] += localNorm2(s.getP()[mu]); 
	//ke_per_site[mu] -= real(trace(s.getP()[mu]*s.getP()[mu]));
      }

      /* Sum up the differences */
      Double KE=zero;
      for(int mu=0; mu < Nd; mu++) { 
	KE += sum(ke_per_site[mu]);
      }

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "mesKE");
      write(xml_out, "KE", sum(KE));
      pop(xml_out);  // pop(mesKE);


      return KE;

      END_CODE();
    }

    //! The Potential Energy 
    Double  mesPE(const AbsFieldState< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >& s) const 
    {
      START_CODE();

      // Self Encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "mesPE");
      // Cycle through all the monomials and compute their contribution
      int num_terms = monomials.size();

      write(xml_out, "num_terms", num_terms);
      Double PE=zero;

      // Caller writes elem rule
      push(xml_out, "PEByMonomials");
      for(int i=0; i < num_terms; i++) 
      {
	push(xml_out, "elem");
	Double tmp;
	tmp=monomials[i]->S(s);
	PE += tmp;
	pop(xml_out); // elem
      }
      pop(xml_out); // PEByMonomials
      pop(xml_out); // pop(mesPE);
      
      END_CODE();
      return PE;
    }


     //------------------------------------------------------------------
     /*
      Driver subroutine to compute all the 3nd and 5th order 
      Poisson Brackets for a given list of monomials

      Paulo Silva, in collaboration with Balint Joo, September 21st, 2010
     */
     //------------------------------------------------------------------
     void mesPB(const AbsFieldState< multi1d<LatticeColorMatrix>,
              multi1d<LatticeColorMatrix> >& s, string& state) const {

       START_CODE();


       typedef multi1d<LatticeColorMatrix>  P;
       typedef multi1d<LatticeColorMatrix>  Q;
 
       QDPIO::cout << "PSILVA: computing PB " << state << endl;

       int nmon = monomials.size();

       QDPIO::cout << "PSILVA: number of monomials = " << nmon << endl;


       multi1d< Handle< TowerArray<LatticeColorMatrix> > > F(nmon);
       multi2d< Handle< TowerArray<LatticeColorMatrix> > > G(nmon, nmon);
       // multi1d< Handle< TowerArray<typename PQTraits<Q>::Base_t> > > F(nmon);
       //multi2d< Handle< TowerArray<typename PQTraits<Q>::Base_t> > > G(nmon, nmon);
       for (int i=0; i< nmon; i++){
	 F[i] = new TowerArray<LatticeColorMatrix>(4);
         
         for(int mu=0; mu < Nd; mu++) { 
	   (*(F[i]))[mu] = zero;
	 }
        
         for (int j=0; j< nmon; j++){
           G[i][j] = new TowerArray<LatticeColorMatrix>(2);

           for(int mu=0; mu < Nd; mu++) { 
            (*G[i][j])[mu] = zero;
	   }
         }
       }

 
      // Computing F towers
      for (int i=0; i< nmon; i++) {

	QDPIO::cout << " PSILVA: Computing F towers: " << i  << endl;
        monomials[i]->dsdq(*(F[i]),s);


        LatticeColorMatrix tet;

        tet=(*(F[i]))[3][0];

        Double rrr;
        rrr=sum(real(trace(tet*adj(tet))));

	QDPIO::cout << "TESTPSILVA: 0  " << rrr << endl;

        tet=(*(F[i]))[3][1];
        rrr=sum(real(trace(tet*adj(tet))));
	QDPIO::cout << "TESTPSILVA: 1  " << rrr << endl;

        tet=(*(F[i]))[3][2];
        rrr=sum(real(trace(tet*adj(tet))));
	QDPIO::cout << "TESTPSILVA: 2  " << rrr << endl;

        tet=(*(F[i]))[3][3];
        rrr=sum(real(trace(tet*adj(tet))));
	QDPIO::cout << "TESTPSILVA: 3  " << rrr << endl;

      }

      //Computing G towers
      for (int i=0; i< nmon; i++) {

        Handle< AbsFieldState<P,Q> > s_new(s.clone());

        for(int mu=0; mu < Nd; mu++) { 
          s_new->getP()[mu] = (*(F[i]))[mu][0];
        }

        for (int j=0; j< nmon; j++){
          
          QDPIO::cout << " PSILVA: Computing G towers: " << i << j  << endl;
	  monomials[j]->dsdq(*(G[i][j]),(*s_new)); 
        }
      }

      //Getting momenta
      P pp;
      pp.resize(Nd);
      pp=s.getP();

       QDPIO::cout << " PSILVA: Got momenta, now computing PB's "  << endl;

       //***** Computation of the Poisson Brackets
       //Declaring...
       multi1d<Double> tst, tttst;
       multi2d<Double> sst, ttsst, sttst;
       multi3d<Double> stsst;

       //Resizing...
       tst.resize(nmon); tttst.resize(nmon); 
       sst.resize(nmon, nmon); 
       ttsst.resize(nmon, nmon);
       sttst.resize(nmon, nmon);
       stsst.resize(nmon, nmon, nmon);

       //Initializing...
       for (int i=0; i<nmon; i++) {
         tst[i]=zero; tttst[i]=zero;
         for (int j=0; j<nmon; j++) {           
           sst[i][j]=zero; ttsst[i][j]=zero; sttst[i][j]=zero;
           for (int k=0; k<nmon; k++) {
             stsst[i][j][k]=zero;
	}}}


       //NEED FOR AN ADDITIONAL FACTOR OF ( -2 ) ... NOT YET IMPLEMENTED

      for(int mu=0; mu < Nd; mu++) {

	 for (int i=0; i<nmon; i++) {

           QDPIO::cout << " PSILVA: computing PBs i = " << i  << endl;

	   tst[i] += Double(-1.0)*sum(real(trace((*F[i])[mu][1]*pp[mu]))); 
           tttst[i] += Double(-1.0)*sum(real(trace((*F[i])[mu][3]*pp[mu])));

	   for (int j=0; j<nmon; j++) {
 
             QDPIO::cout << " PSILVA: computing PBs j = " << j  << endl;

             LatticeColorMatrix com_i, com_j, com1p, com2p;
           
	     sst[i][j] += Double(1.0)*sum(real(trace((*F[i])[mu][0]*(*F[j])[mu][0])));

             ttsst[i][j] += Double(2.0)*sum(real(trace((*F[i])[mu][1]*(*F[j])[mu][1])));
             ttsst[i][j] += Double(1.0)*sum(real(trace((*F[i])[mu][0]*(*F[j])[mu][2])));
             ttsst[i][j] += Double(1.0)*sum(real(trace((*F[i])[mu][2]*(*F[j])[mu][0]))); 

 
	     com1p = (*F[i])[mu][0]*(*F[j])[mu][1]-(*F[j])[mu][1]*(*F[i])[mu][0];
	     com2p = (*F[i])[mu][1]*(*F[j])[mu][0]-(*F[j])[mu][0]*(*F[i])[mu][1];  

             sttst[i][j] += Double(-2.0)*sum(real(trace(com1p*pp[mu])));
             sttst[i][j] += Double(1.0)*sum(real(trace(com2p*pp[mu])));

	     com_i=(*F[i])[mu][0]*pp[mu]-pp[mu]*(*F[i])[mu][0];
             com_j=(*F[j])[mu][0]*pp[mu]-pp[mu]*(*F[j])[mu][0];

 	     sttst[i][j] += Double(-1.0)*sum(real(trace(com_i*com_j))); 
             sttst[i][j] += Double(1.0)*sum(real(trace((*F[i])[mu][0]*(*F[j])[mu][2])));
 	     sttst[i][j] += Double(-2.0)*sum(real(trace((*F[i])[mu][1]*(*F[j])[mu][1])));

 	     for (int k=0; k<nmon; k++) {

                QDPIO::cout << " PSILVA: computing PBs k = " << k  << endl;

                LatticeColorMatrix aux1, aux2, aux3, aux4;
                aux1=(*F[j])[mu][0];
                aux2=(*G[i][k])[mu][1];

                aux3=(*F[k])[mu][0];
                aux4=(*G[i][j])[mu][1];

                stsst[i][j][k] += Double(-1.0)*sum(real(trace(aux1*aux2+aux3*aux4)));   

 	     }
	   }
	 }


       }//mu



      //NOW, OUTPUT THIS HUGE AMOUNT OF STUFF...
       for (int i=0; i<nmon; i++) {    
         QDPIO::cout << " PPGFsep " << state << " tst" << i << "  " << tst[i] << endl;
         QDPIO::cout << " PPGFsep " << state << " tttst" << i << "  " << tttst[i] << endl;
         for (int j=0; j<nmon; j++) {           
           QDPIO::cout << " PPGFsep " << state << " sst" << i << j << "  " << sst[i][j] << endl;
	    QDPIO::cout << " PPGFsep " << state << " ttsst" << i << j << "  " << ttsst[i][j] << endl;
           QDPIO::cout << " PPGFsep " << state << " sttst" << i << j << "  " << sttst[i][j] << endl;
           for (int k=0; k<nmon; k++) {
             QDPIO::cout << " PPGFsep " << state << " stsst" << i << j << k << "  " << stsst[i][j][k] << endl;
	}}} 





      END_CODE();

 
     }



  private:
    //! Convenience 
    typedef ExactMonomial< multi1d<LatticeColorMatrix>, 
			           multi1d<LatticeColorMatrix> >  ExactMon;

    //! This creates the hamiltonian. It is similar to the 
    void create(const multi1d<std::string>& monomial_ids);
    

    multi1d< Handle<ExactMon> >  monomials;

    
  };

}

#endif
