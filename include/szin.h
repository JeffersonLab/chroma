// $Id: szin.h,v 1.1 2002-12-16 06:56:49 edwards Exp $
//
// Include file for test suite

#include <qdp.h>
#include <geom.h>

using namespace QDP;

enum Reunitarize {REUNITARIZE, REUNITARIZE_ERROR, REUNITARIZE_LABEL};
enum Exponentiate {TWELWTH_ORDER, EXACT};
enum Sources {POINT_SOURCE, WALL_SOURCE, SHELL_SOURCE, BNDST_SOURCE, POINT_AND_BNDST_SOURCE, 
	      SHELL_AND_BNDST_SOURCE, POINT_AND_SHELL_AND_BNDST_SOURCE};
enum Sinks {POINT_SINK, WALL_SINK, POINT_AND_WALL_SINK, SHELL_SINK, POINT_AND_SHELL_SINK, 
	    BNDST_SINK, POINT_AND_BNDST_SINK, SHELL_AND_BNDST_SINK, POINT_AND_SHELL_AND_BNDST_SINK};



void junk(LatticeGauge& b3, const LatticeGauge& b1, const LatticeGauge& b2, const Subset& s);
void MesPlq(const multi1d<LatticeGauge>& u, Double& w_plaq, Double& s_plaq, 
	    Double& t_plaq, Double& link);
void mesons(const LatticePropagator& quark_prop_1, const LatticePropagator& quark_prop_2, 
	    multi2d<Real>& meson_propagator, 
	    int num_mom, const multi1d<int>& t_source, 
	    int j_decay);
void baryon(LatticePropagator& quark_propagator, 
	    multi2d<Complex>& barprop, 
	    const multi1d<int>& t_source, int j_decay, int bc_spec);
void dslash(LatticeFermion& chi, const multi1d<LatticeColorMatrix>& u, const LatticeFermion& psi, int isign, int cb);

void mciter(multi1d<LatticeGauge>& u, int n_over, int nheat, int &ntrials, int &nfails);

void reunit(LatticeColorMatrix& a);
void reunit(LatticeColorMatrix& xa, LatticeBoolean& bad, enum Reunitarize ruflag, int& numbad);
void su3over(LatticeColorMatrix& u, const LatticeColorMatrix& w, int su2_index);
void su3hb(LatticeColorMatrix& u, const LatticeColorMatrix& w, 
	   int su2_index, int nheat, int& ntrials, int& nfails);

void taproj(LatticeColorMatrix& ds_u);
void eeu1(LatticeColorMatrix& a);
void eesu2(LatticeColorMatrix &a);
void expsu3(LatticeColorMatrix &a, enum Reunitarize opt);
void expm12(LatticeColorMatrix& a);
void expmat(LatticeColorMatrix& a, enum Exponentiate opt);
void srcfil(LatticeFermion& a, const multi1d<int>& coord, int colour_index, int spin_index);
void walfil(LatticeFermion& a, int slice, int mu, int colour_index, int spin_index);

void HotSt(multi1d<LatticeColorMatrix>& u);
void rgauge(multi1d<LatticeColorMatrix>& u, LatticeColorMatrix& g);


//void CREATE_NAMELIST(const char *s, ...);
//inline void WRITE_NAMELIST(...) {}
#define CREATE_NAMELIST(a,b) 

// #define WRITE_NAMELIST(grp,...) Write(cout, #grp , ##__VA_ARGS__)


const float fuzz = 1.0e-5;
const float twopi = 6.283185307179586476925286;

#define START_CODE(a)
#define END_CODE()

#define TO_REAL(a) float(a)


#include <defs.h>
#include "common_declarations.h"
