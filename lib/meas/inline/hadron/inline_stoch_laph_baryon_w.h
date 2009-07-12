// -*- C++ -*-
// $Id: inline_stoch_laph_baryon_w.h,v 3.2 2009-07-12 00:44:43 jbulava Exp $
/*! \file
 * \brief Inline measurement of stochastic source and sink functions 
 * for baryons
 */

#ifndef __inline_stoch_laph_baryon_h__
#define __inline_stoch_laph_baryon_h__

#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
	/*! \ingroup inlinehadron */
	namespace InlineStochLapHBaryonEnv 
	{
		extern const std::string name;
		bool registerAll();


		//! Parameter structure
		/*! \ingroup inlinehadron */
		struct Params
		{
			Params();
			Params(XMLReader& xml_in, const std::string& path);
			void writeXML(XMLWriter& xml_out, const std::string& path);

			unsigned long      frequency;

			//Noise ids
			//T0's 
			//cfg info
			//dilution scheme
			//smearing info (link and quark)

			struct Param_t
			{

				multi1d<BaryonOperator> bops;               /*!< Baryon operator 
																											information 
																											(see baryon_operator.h)*/

				GroupXML_t          link_smearing;         /*!< link smearing xml */

			};

			struct NamedObject_t
			{
		
				//Names of files that were output of STOCH_LAPH_QUARK
			
			};

			Param_t        param;      /*!< Parameters */    
			NamedObject_t  named_obj;  /*!< Named objects */
		
		};


		//! Inline measurement of stochastic group baryon operators
		/*! \ingroup inlinehadron */
		class InlineMeas : public AbsInlineMeasurement 
		{
			public:
				~InlineMeas() {}
				InlineMeas(const Params& p) : params(p) {}
				InlineMeas(const InlineMeas& p) : params(p.params) {}

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


		struct Sink_qqq_t
		{
			struct DilutionComponent_t
			{
				multi1d<DComplex> time;   //Will be length Lt	
			
				DilutionComponent_t& operator*(const DComplex& coeff)
				{
					DilutionComponent_t result;

					result.time = time * coeff;
				
					return result;
				}
			
			};
		
		
			multi3d<DilutionComponent_t> dilutions;
		
		};
		/*Accessed as follows 
		Sink_qqq_t snk;
		snk.dilutions(d1,d2,d3).time(t);
		*/

		struct Source_qqq_t
		{
			Source_qqq_t(int dil_size) : dilutions.resize(dil_size, dil_size, 
					dil_size) {}

			multi3d<DComplex> dilutions;
		};
		/*Acessed as follows:
		Source_t src;
		src.dilutions(d1,d2,d3);
		*/


		//vector of momenta:   mutli1d< multi1d<int> > momenta( 
		//array of results: multi1d<BaronOpSourceSink_t> bresults( 
		//BaryonOpSourceSink_t(dil_size, nt), momenta.size() ) 
		 

		struct BaryonOpSourceSink_t
		{
			
			BaryonOpSourceSink_t(int dil_size, int nt) : src(dil_size), 
			snk(dil_size, nt)
			{
			}
			
			Source_qqq_t src;

			Sink_qqq_t snk;
		
			BaryonOpSourceSink_t& operator*(const DComplex& coeff)
			{

				BaryonOpSourceSink_t result;

				//Do src first
				result.src.dilutions = src.dilutions * coeff;

				//Sink
				result.snk.dilutions = snk.dilutions * coeff;

				return result;
			}

			BaryonOpSourceSink_t& operator+=(const BaryonOpSourceSink_t& rhs)
			{

				//Check Sizes
				if 
			}
		
		};

		


		
	} // namespace InlineStochGroupBaryonEnv 
}

#endif
