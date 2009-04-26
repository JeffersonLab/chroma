#ifndef SCALAR_LOCAL_S_H
#define SCALAR_LOCAL_S_H

#include "meas/hadron/hadron_corr_s.h"
#include "meas/hadron/stag_propShift_s.h"

namespace Chroma {

  class staggered_hadron_corr;

  class staggered_local_scalar : public staggered_hadron_corr {

  public:
    void compute(LatticeStaggeredPropagator& quark_prop_A,
		 LatticeStaggeredPropagator& quark_prop_B, int j_decay);

    void compute(multi1d<LatticeStaggeredPropagator>& quark_props, int j_decay) {
      if (quark_props.size() != 2) {
	QDPIO::cerr << "staggered_local_scalar::compute(multi1d<LatticeStaggeredPropagator>&, int) : Array arg must contain 2 props\n";
	QDP_abort(1);
      } else {
	compute(quark_props[0], quark_props[1], j_decay);
      }
    }

    staggered_local_scalar(int t_len, const multi1d<LatticeColorMatrix>& u_in,
			   Stag_shift_option type_of_shift_in = SYM_GAUGE_INVAR)
      : staggered_hadron_corr(t_len, no_pions, u_in, type_of_shift_in)
      {
	outer_tag = "Scalar";
	inner_tag = "Sc";

	tag_names.resize(no_pions);
	for (int i=0;i<no_pions;++i) {
	  ostringstream tag;
	  tag << "re_sc" << i;
	  tag_names[i] = tag.str();
	}
      }

    virtual ~staggered_local_scalar() {};

  protected:
    
  private:
    static const int no_pions = 1;

  };

} //end namespace

#endif
