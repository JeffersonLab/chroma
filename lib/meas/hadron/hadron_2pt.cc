// $Id: hadron_2pt.cc,v 1.1 2007-05-09 17:19:44 edwards Exp $
/*! \file
 *  \brief Construct hadron correlators
 */

#include "meas/hadron/hadron_2pt.h"
#include "meas/inline/io/named_objmap.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{

  // Get source location
  // Default version
  multi1d<int>
  HadronCorrelator::getTSrce(const multi1d<ForwardProp_t>& forward_headers) const
  {
    multi1d<int> t_srce = forward_headers[0].source_header.getTSrce();
    int j_decay = forward_headers[0].source_header.j_decay;

    for(int loop=1; loop < forward_headers.size(); ++loop)
    {
      multi1d<int> t_srce_b = forward_headers[loop].source_header.getTSrce();

      // Bummer, I wish qdp++ had a multi1d.operator!=()
      bool same = true;
      for(int i=0; i < t_srce.size(); ++i)
      {
	if (t_srce_b[i] != t_srce[i]) 
	  same = false;
      }
      
      if (! same)
      {
	QDPIO::cerr << __func__ << ": the t_srce in the forward props are not all equal"
		    << endl;
	QDP_abort(1);
      }

      int j_decay_b = forward_headers[loop].source_header.j_decay;

      if (j_decay != j_decay_b)
      {
	QDPIO::cerr << __func__ << ": the decay_dir in the forward props are not all equal"
		    << endl;
	QDP_abort(1);
      }
    }

    return t_srce;
  }


  //
  // Extract quark prop headers
  // Default version
  //
  ForwardProp_t
  HadronCorrelator::readPropHeader(const std::string& prop_id) const
  {
    ForwardProp_t header;

    try
    {
// No snarfing 
//	ret.prop = TheNamedObjMap::Instance().getData<T>(prop_id);
	
      // Snarf the source info. This is will throw if the source_id is not there
      XMLReader prop_file_xml, prop_record_xml;
      TheNamedObjMap::Instance().get(prop_id).getFileXML(prop_file_xml);
      TheNamedObjMap::Instance().get(prop_id).getRecordXML(prop_record_xml);
   
      // Try to invert this record XML into a ChromaProp struct
      read(prop_record_xml, "/Propagator", header);
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << __func__ << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << __func__ << ": caught error = " << e 
		  << endl;
      QDP_abort(1);
    }

    return header;
  }


  //! Project onto fixed momenta
  std::list<Hadron2PtCorrs_t> HadronCorrelator::project(const std::list<Hadron2PtContract_t>& had_list, 
							const Hadron2PtCorrParams_t& params) const
  {
    SftMom phases(params.mom2_max,
		  params.t_srce,
		  params.mom_origin,
		  params.avg_equiv_mom,
		  params.decay_dir);

    int length = phases.numSubsets();

    std::list<Hadron2PtCorrs_t> corrs;

    // Run over the input list, 
    for(std::list<Hadron2PtContract_t>::const_iterator had_ptr= had_list.begin(); 
	had_ptr != had_list.end(); 
	++had_ptr)
    {
      const Hadron2PtContract_t& had_cont = *had_ptr;
      multi2d<DComplex> hsum(phases.sft(had_cont.corr));   // slow fourier-transform

      // Copy onto output structure
      Hadron2PtCorrs_t  had_corrs;
      had_corrs.xml = had_cont.xml;

      for(int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) 
      {
	Hadron2PtCorrs_t::Mom_t  had_mom;
	had_mom.mom = phases.numToMom(sink_mom_num);
     
	had_mom.corr.resize(length);   // QUESTION: DO WE WANT TO CHANGE THE ORIGIN CONVENTION?? YES!!
	for (int t=0; t < length; ++t) 
	{
//        int t_eff = (t - t0 + length) % length;
	  had_mom.corr[t] = hsum[sink_mom_num][t];
	}
	
	had_corrs.corrs.push_back(had_mom);
      }

      corrs.push_back(had_corrs);
    }

    return corrs;
  }


} // end namespace Chroma
