// $Id: hadron_2pt.cc,v 1.4 2007-06-12 16:09:37 edwards Exp $
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
    Handle<HadronContractResult_t> output(new HadronContractResult_t);
    output->xml << xml;
    output->xml_regres << xml_regres;

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
  Hadron2PtCorr::project(const std::list< Handle<Hadron2PtContract_t> >& had_list, 
			 const SftMomParams_t& params) const
  {
    std::list< Handle<HadronContractResult_t> > corrs;

    SftMom phases(params);

    int length = phases.numSubsets();

    // Run over the input list, 
    for(std::list< Handle<Hadron2PtContract_t> >::const_iterator had_ptr= had_list.begin(); 
	had_ptr != had_list.end(); 
	++had_ptr)
    {
      const Hadron2PtContract_t& had_cont = **had_ptr;
      multi2d<DComplex> hsum(phases.sft(had_cont.corr));   // slow fourier-transform

      // Copy onto output structure
      Hadron2PtCorrs_t  had_corrs;
      had_corrs.xml << had_cont.xml;

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

      // Use the zero momentum data as a regression output
      // NOTE: the name of this group is not important, and you can jam whatever
      // you want into here. It is used for the regressions to latch onto something
      // from the output since it is all in binary
      push(had_corrs.xml_regres, "Hadron2Pt");
      write(had_corrs.xml_regres, "ZeroMom", hsum[0]);
      pop(had_corrs.xml_regres);

      // Serialize the object and put it onto the end of the list
      corrs.push_back(had_corrs.serialize());
    }

    return corrs;
  }


} // end namespace Chroma
