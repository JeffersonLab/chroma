// $Id: hadron_2pt.cc,v 1.3 2007-06-10 14:49:06 edwards Exp $
/*! \file
 *  \brief Construct hadron 2pt correlators
 */

#include "meas/hadron/hadron_2pt.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{

  // Serialize the structure
  Handle<HadronContractResult_t> Hadron2PtCorrs_t::serialize() const
  {
    Handle<HadronContractResult_t> output;
    output->xml = xml;

    // Run over the momentum list 
    for(std::list<Hadron2PtCorrs_t::Mom_t>::const_iterator mom_ptr= corrs.begin(); 
	mom_ptr != corrs.end(); 
	++mom_ptr)
    {
      write(output->bin, mom_ptr->mom);
      write(output->bin, mom_ptr->corr);
    }

    return output;
  }


  //! Project onto fixed momenta
  std::list< Handle<HadronContractResult_t> > 
  Hadron2PtCorr::project(const std::list<Hadron2PtContract_t>& had_list, 
			 const Hadron2PtCorrParams_t& params) const
  {
    std::list< Handle<HadronContractResult_t> > corrs;

    SftMom phases(params.mom2_max,
		  params.t_srce,
		  params.mom_origin,
		  params.avg_equiv_mom,
		  params.decay_dir);

    int length = phases.numSubsets();

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

      corrs.push_back(had_corrs.serialize());
    }

    return corrs;
  }


} // end namespace Chroma
