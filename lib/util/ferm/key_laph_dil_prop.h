// -*- C++ -*-
// $Id: key_laph_dil_prop.h,v 1.2 2009-07-09 02:13:21 jbulava Exp $
/*! \file
 * \brief Key for propagator colorvector matrix elements
 */

#ifndef __key_laph_dil_prop_h__
#define __key_laph_dil_prop_h__

#include "chromabase.h"
#include "util/ferm/key_val_db.h"

namespace Chroma
{
  //----------------------------------------------------------------------------
  /*!
   * \ingroup ferm
   * @{
   */
  //! Prop operator

	struct KeyLaphDilutedProp_t
		{
			std::string src_or_snk; //Is this a source or solution? Will be either
															//"SRC" or  "SNK"

			int t0; //If the above string is 'SRC', this will always be 0. If it is 
							//"SNK", this will be the t0 on which it was created

			int spin_dil; //The spin dilution component

			int evec_dil; //The eigenvector dilution component

			int noise_src; //The noise index

		};

		bool operator<(const KeyLaphDilutedProp_t& a, 
				const KeyLaphDilutedProp_t& b);

		void write(BinaryWriter& bin, const KeyLaphDilutedProp_t& param);
		
		void read(BinaryReader& bin, KeyLaphDilutedProp_t& param);

		void read(XMLReader& xml, const std::string& path, 
			KeyLaphDilutedProp_t& param);

		void write(XMLWriter& xml, const std::string& path, 
			const KeyLaphDilutedProp_t& param);

 
} // namespace Chroma

#endif
