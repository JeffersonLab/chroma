// -*- C++ -*-
// $Id: key_laph_dil_prop.cc,v 1.2 2009-07-09 02:13:21 jbulava Exp $
/*! \file
 * \brief Key for propagator colorvector matrix elements
 */

#include "util/ferm/key_laph_dil_prop.h"

namespace Chroma 
{ 
 
		bool operator<(const KeyLaphDilutedProp_t& a, 
				const KeyLaphDilutedProp_t& b)
		{
			
			bool ret; 
			
			multi1d<int> lgaa(3);
			lgaa[0] = a.spin_dil;
			lgaa[1] = a.evec_dil;
			lgaa[2] = a.noise_src;
			//lgaa[3] = a.src_or_snk;

			multi1d<int> lgbb(3);
			lgbb[0] = b.spin_dil;
			lgbb[1] = b.evec_dil;
			lgbb[2] = b.noise_src;
			//lgbb[3] = b.src_or_snk;
   

			if (lgaa == lgbb)
				ret = (a.src_or_snk < b.src_or_snk);
			else 
				ret = (lgaa < lgbb);

			return ret;
		}

		void write(BinaryWriter& bin, const KeyLaphDilutedProp_t& param)
		{
			write(bin, param.src_or_snk);
			write(bin, param.spin_dil);
			write(bin, param.evec_dil);
			write(bin, param.noise_src);
		}

		void read(BinaryReader& bin, KeyLaphDilutedProp_t& param)
		{
			read(bin, param.src_or_snk, 32);
			read(bin, param.spin_dil);
			read(bin, param.evec_dil);
			read(bin, param.noise_src);
		}

		//! KeyPropColorVec reader
  void read(XMLReader& xml, const std::string& path, 
			KeyLaphDilutedProp_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "src_or_snk", param.src_or_snk);
    read(paramtop, "spin_dil", param.spin_dil);
    read(paramtop, "evec_dil", param.evec_dil);
    read(paramtop, "noise_src", param.noise_src);
  
	}

  // KeyPropColorVec writer
  void write(XMLWriter& xml, const std::string& path, 
			const KeyLaphDilutedProp_t& param)
  {
    push(xml, path);

    write(xml, "src_or_snk", param.src_or_snk);
    write(xml, "spin_dil", param.spin_dil);
    write(xml, "evec_dil", param.evec_dil);
    write(xml, "noise_src", param.noise_src);

    pop(xml);
  }



} // namespace Chroma
