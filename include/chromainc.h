// $Id: chromainc.h,v 1.1 2003-02-16 04:02:22 edwards Exp $
//
// Include file that includes all the include files.

#include "qdp.h"

#include "common_declarations.h"
#include "primitives.h"

#include "actions/actions.h"
#include "meas/meas.h"
#include "update/update.h"
#include "util/util.h"

using namespace QDP;

enum Exponentiate {TWELWTH_ORDER, EXACT};
enum Sources {POINT_SOURCE, WALL_SOURCE, SHELL_SOURCE, BNDST_SOURCE, POINT_AND_BNDST_SOURCE, 
	      SHELL_AND_BNDST_SOURCE, POINT_AND_SHELL_AND_BNDST_SOURCE};
enum Sinks {POINT_SINK, WALL_SINK, POINT_AND_WALL_SINK, SHELL_SINK, POINT_AND_SHELL_SINK, 
	    BNDST_SINK, POINT_AND_BNDST_SINK, SHELL_AND_BNDST_SINK, POINT_AND_SHELL_AND_BNDST_SINK};



void mesons(const LatticePropagator& quark_prop_1, const LatticePropagator& quark_prop_2, 
	    multi2d<Real>& meson_propagator, 
	    int num_mom, const multi1d<int>& t_source, 
	    int j_decay);
void baryon(LatticePropagator& quark_propagator, 
	    multi2d<Complex>& barprop, 
	    const multi1d<int>& t_source, int j_decay, int bc_spec);
void dslash(LatticeFermion& chi, const multi1d<LatticeColorMatrix>& u, const LatticeFermion& psi, int isign, int cb);

void mciter(multi1d<LatticeColorMatrix>& u, int n_over, int nheat, int &ntrials, int &nfails);

void su3over(LatticeColorMatrix& u, const LatticeColorMatrix& w, int su2_index);
void su3hb(LatticeColorMatrix& u, const LatticeColorMatrix& w, 
	   int su2_index, int nheat, int& ntrials, int& nfails);

void walfil(LatticeFermion& a, int slice, int mu, int colour_index, int spin_index);

void eeu1(LatticeColorMatrix& a);
void eesu2(LatticeColorMatrix &a);
void expsu3(LatticeColorMatrix &a, enum Reunitarize opt);
void expmat(LatticeColorMatrix& a, enum Exponentiate opt);

void rgauge(multi1d<LatticeColorMatrix>& u, LatticeColorMatrix& g);

