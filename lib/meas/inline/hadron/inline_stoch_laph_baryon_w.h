// -*- C++ -*-
// $Id: inline_stoch_laph_baryon_w.h,v 3.3 2009-07-13 03:59:49 jbulava Exp $
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

				int nnoises;         //The number of noise source to be used (must be 
														//>= 3 

			};

			struct NamedObject_t
			{
		
				//Names of files that were output of STOCH_LAPH_QUARK
		
				multi1d<std::string> quark_files; //The db files where the diluted
																					//sources and sinks are stored
																					
				std::string baryon_file; //The output db file that stores the 
																 //Operator source and sink functions.
				
				std::string gauge_id;
			
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
			
			Sink_qqq_t(int dil_size, int nt) : dilutions.resize(dil_size, 
					dil_size, dil_size) 
			{
				

				for (int d1 = 0 ; d1 < dil_size ; ++d1)
					for (int d2 = 0 ; d2 < dil_size ; ++d2)
						for (int d3 = 0 ; d3 < dil_size ; ++d3)
						{
							dilutions(d1, d2, d3).time.resize(nt);

							for(int t = 0 ; t < nt ; ++t)
								dilutions(d1, d2, d3).time[t] = Complex(0.0);
				
						}

			}	
			
			
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
					dil_size) 
			{
				for (int d1 = 0 ; d1 < dil_size ; ++d1)
					for (int d2 = 0 ; d2 < dil_size ; ++d2)
						for (int d3 = 0 ; d3 < dil_size ; ++d3)
						{
							dilutions(d1, d2, d3) = Complex(0.0);
				
						}
			}

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
			snk(dil_size, nt) {}
			
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

				src.dilutions += rhs.src.dilutions;

				snk.dilutions += rhs.snk.dilutions;
				 
				return this;
			}
		
		};

	
		struct KeyBaryonOpSourceSink_t
		{
			int op_num;
			
			multi1d<int> noises; 
			int t0; 

		};

			bool operator<(const KeyBaryonOpSourceSink_t& a, 
				const KeyBaryonOpSourceSink_t& b);


			void write(BinaryWriter& bin, const KeyBaryonOpSourceSink_t& param);
		
			void read(BinaryReader& bin, KeyBaryonOpSourceSink_t& param);

			void read(XMLReader& xml, const std::string& path, 
				KeyBaryonOpSourceSink_t& param);

		void write(XMLWriter& xml, const std::string& path, 
			const KeyBaryonOpSourceSink_t& param);
		
	} // namespace InlineStochGroupBaryonEnv 
}

#endif
