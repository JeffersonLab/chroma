// -*- C++ -*-
/*! \file
 *
 * \brief 
 *
 * Helper structures that are used for operators inserted in genprops and NPR
 * 
 */

#ifndef  __meson_ops_h__
#define  __meson_ops_h__

// Include everything...
#include "chroma.h"


namespace Chroma 
{ 
  namespace MesonOps{
    struct Operators_t{
      std::string name;
      int gamma ;
    };

    
    struct State_t{
      std::string name;
      multi1d<std::string> flavor ;
      std::string db ;
      std::string wavefunc_file ;
    };

    
    void read(XMLReader& xml, const std::string& path, Operators_t& op);
    void write(XMLWriter& xml, const std::string& path, const Operators_t& op);

    void read(XMLReader& xml,const std::string& path, State_t& s);
    void write(XMLWriter& xml, const std::string& path, const State_t& s);
    
  }

  namespace MesSpec{
    
    void contract(LatticeComplex& latC, 
		  const LatticePropagator& aq,
		  const LatticePropagator& q,
		  int snk,
		  int src);
    
    
    class Key{
    public:
      int P    ; // Parity
      int J    ; // 2 x Total spin (only spin 1/2 and spin 3/2 will be used... )
      int Jz   ; // 2 x Z component of spin
      int I    ; // 2 x total isospin
      int I3   ; // 2 x 3rd component of isospin
      int S    ; // strangness 
      int C    ; // charm (being generic here....
      int B    ; // bottom (being more generic here....
      
      Key():P(1),J(0),Jz(0),I(0),I3(0),S(0),C(0),B(0){}// default constructor
      Key(int p_, int j_, int jz_,
	  int i_, int i3_,int s_,int c_, int b_, int id_):
	P(p_),J(j_),Jz(jz_),I(i_),I3(i3_),S(s_),C(c_),B(b_){}

      // NOTE THAT THE KeyQQQBlock serialization in MULTIBARSPEC
      // HAS I Iz before J and Jz
      multi1d<int> serialize(){
	multi1d<int> ret(8);
	ret[0]=P; 
	ret[1]=J; ret[2]=Jz;
	ret[3]=I; ret[4]=I3; 
	ret[5]=S; 
	ret[6]=C;
	ret[7]=B;
	return ret ;
      } ;
      
      ~Key(){}
    };
    void read(XMLReader& xml, const std::string& path, Key& k) ;
    void write(XMLWriter& xml, const std::string& path, const Key& k);
    
  }// namespace MesSpec



}
#endif
