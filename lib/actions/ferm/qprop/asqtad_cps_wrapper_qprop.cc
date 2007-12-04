#include "chromabase.h"

#include "actions/ferm/qprop/asqtad_cps_wrapper_qprop.h"
#include "actions/ferm/invert/syssolver_cg_params.h"
#include "util/gauge/stag_phases_s.h"


extern "C" {

#define QLA_Precision 'F'
#define QLA_Colors 3

#include <qdp.h>
#include <qdp_types.h>
#include <qdp_common.h>
#include <qdp_f3_generic.h>
#include <qdp_f3.h>
#include <qop_qdp.h>
#include <qmp.h>

#include <qla.h>
#include <qla_complex.h>


}




/**************************************/
/** level 3 interface code  **/


//
//  utility routines
//

//
//  from test_common.c
//

QLA_Real
get_plaq(QDP_ColorMatrix *link[])
{
  int mu, nu;
  QLA_Real plaq;
  QDP_ColorMatrix *temp1, *temp2, *temp3, *temp4;

#ifdef LOCAL_SUM
  QDP_Real *treal1, *treal2;
  treal1 = QDP_create_R();
  treal2 = QDP_create_R();
  QDP_R_eq_zero(treal2, QDP_all);
#else
  QLA_Real tplaq;
  plaq = 0;
#endif

  temp1 = QDP_create_M();
  temp2 = QDP_create_M();
  temp3 = QDP_create_M();
  temp4 = QDP_create_M();

  // http://usqcd.jlab.org/usqcd-docs/qdp/qdpc.html/Functions-involving-shifts.html#Functions-involving-shifts

  for(mu=0; mu<QDP_ndim()-1; ++mu) {
    for(nu=mu+1; nu<QDP_ndim(); ++nu) {

      QDP_M_eq_sM(temp1, link[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
      QDP_M_eq_sM(temp2, link[mu], QDP_neighbor[nu], QDP_forward, QDP_all);

      QDP_M_eq_Ma_times_M(temp3, link[nu], link[mu], QDP_all);

      QDP_M_eq_M_times_M(temp4, temp3, temp1, QDP_all);
      QDP_discard_M(temp1);

#ifdef LOCAL_SUM
      QDP_R_eq_re_M_dot_M(treal1, temp2, temp4, QDP_all);
      QDP_discard_M(temp2);
      QDP_R_peq_R(treal2, treal1, QDP_all);
#else
      QDP_r_eq_re_M_dot_M(&tplaq, temp2, temp4, QDP_all);
      QDP_discard_M(temp2);
      plaq += tplaq;
#endif

    }
  }

#ifdef LOCAL_SUM
  QDP_r_eq_sum_R(&plaq, treal2, QDP_all);
  QDP_destroy_R(treal1);
  QDP_destroy_R(treal2);
#endif

  QDP_destroy_M(temp1);
  QDP_destroy_M(temp2);
  QDP_destroy_M(temp3);
  QDP_destroy_M(temp4);

  // normalisation of plaquette <p> = 1 for cold config.
  return plaq/(3.0*0.5*QDP_ndim()*(QDP_ndim()-1)*QDP_volume());
}



//
// Convert the gaugefield from chroma/qdp++ format 
// into qdp format.
//
QDP_ColorMatrix ** QDP_create_gauge_from_chroma (
					       multi1d<LatticeColorMatrix> & u_chroma)
{
  QDP_ColorMatrix **u;
  
  const int ndim = 4 ;
  u = (QDP_ColorMatrix **) malloc(ndim*sizeof(QDP_ColorMatrix *));
  for(int i=0; i<ndim; i++) u[i] = QDP_create_M();
  // more work --- check error codes

  QLA_ColorMatrix *tmp;
  tmp = 
    (QLA_ColorMatrix *) malloc(QDP_sites_on_node*sizeof(QLA_ColorMatrix));

  for(int i=0; i<ndim; i++)
    {
      for(int site=0 ; site < QDP_sites_on_node ; ++site)
	{
	  for(int ic = 0 ; ic < 3 ; ++ic) 
	    for(int jc = 0 ; jc < 3 ; ++jc) 
	      {

		Real rrr = u_chroma[i].elem(site).elem().elem(ic,jc).real() ; 
		Real iii = u_chroma[i].elem(site).elem().elem(ic,jc).imag() ; 
		
		QLA_Real zre  = toFloat(rrr) ;
		QLA_Real zim  = toFloat(iii) ; 
		QLA_Complex z ; 
		QLA_C_eq_R_plus_i_R(&z,&zre,&zim) ;

		//		cout << "" << ic << " , " << jc << " = " << zre << " " << zim << endl ;
		QLA_M_eq_elem_C(&tmp[site],&z,ic,jc) ; 
	      }

	  //	  tmp[site]
	}
      // convert to qdp format
      QDP_insert_M(u[i], tmp,QDP_all) ; 
    }

  free(tmp) ; 


  // convert to QOP format.
  QLA_Real plaq;
  plaq = get_plaq(u);
  printf("plaquette = %g\n", plaq);

  return u ;
}

//
// Convert the colorvector from chroma/qdp++ format 
// into qdp format.
//
void convert_chroma_to_qdp(QDP_ColorVector *out,  
const LatticeStaggeredFermion & in ) 
{
  QLA_ColorVector *tmp;
  tmp = 
    (QLA_ColorVector *) malloc(QDP_sites_on_node*sizeof(QLA_ColorVector));

  for(int site=0 ; site < QDP_sites_on_node ; ++site)
    {

      for(int ic = 0 ; ic < 3 ; ++ic) 
	{
	  QLA_Real zre  ; 

	  Real rrr = in.elem(site).elem(0).elem(ic).real() ; 
	  Real iii = in.elem(site).elem(0).elem(ic).imag() ; 
		
	  zre  = toFloat(rrr) ;
	  QLA_Real zim  = toFloat(iii) ; 
	  QLA_Complex z ; 
	  QLA_C_eq_R_plus_i_R(&z,&zre,&zim) ;

	  QLA_V_eq_elem_C(&tmp[site],&z,ic) ; 

	  //	  cout << "convert_chroma_to_qdp-DEBUG " << site << " " << ic << " " << zre << "  " <<  zim << "\n" ;
	}

    } // loop over sites

  // convert to qdp format
  QDP_insert_V(out, tmp,QDP_all) ; 


  free(tmp) ; 

}



//
// Convert the colorvector in qdp into  chroma/qdp++ format 
//
//
void convert_qdp_to_chroma(LatticeStaggeredFermion & out,
			   QDP_ColorVector *in) 
{
  Complex zz ; 
  QLA_Complex z ; 

  QLA_ColorVector *tmp;
  tmp = QDP_expose_V(in) ; 
  // this is done at the pointer level

  for(int site=0 ; site < QDP_sites_on_node ; ++site)
    {
      for(int ic = 0 ; ic < 3 ; ++ic) 
	{
	  QLA_C_eq_elem_V(&z,&tmp[site],ic) ; 

	  QLA_Real zre  ; 
	  zre = QLA_real(z) ;

	  QLA_Real zim  ; 
	  zim = QLA_imag(z) ;

	  Real zre_chroma , zim_chroma ;
	  zre_chroma = (Real)  zre ;
	  zim_chroma = (Real)  zim ;

	  out.elem(site).elem(0).elem(ic).real() = toFloat(zre_chroma) ;
	  out.elem(site).elem(0).elem(ic).imag() = toFloat(zim_chroma) ;


	  //	  cout << "convert_qdp_to_chroma-DEBUG " << site << " " << ic << " " << zre_chroma << " " << zim_chroma   << "\n" ;

	} // loop over color

    } // loop over sites

}


//
// set the coefficients from the MILC code.
// This routine came from 
//           generic_ks/quark_stuff.c
//           generic_ks/imp_actions/asqtad_action.h
//

void load_qop_asqtad_coeffs(QOP_asqtad_coeffs_t *c, QLA_Real weight)
{
  QLA_Real ferm_epsilon;
#define MAX_BASIC_PATHS 6
  static QLA_Real act_path_coeff[MAX_BASIC_PATHS] = {
    ( 1.0/8.0)+(6.0/16.0)+(1.0/8.0),        /* one link */
    /*One link is 1/8 as in fat7 +3/8 for Lepage + 1/8 for Naik */
    (-1.0/24.0),                 /* Naik */
    (-1.0/8.0)*0.5,              /* simple staple */
    ( 1.0/8.0)*0.25*0.5,         /* displace link in two directions */
    (-1.0/8.0)*0.125*(1.0/6.0),  /* displace link in three directions */
    (-1.0/16 ),                  /* Correct O(a^2) errors */
  };


  /* Load path coefficients from table */
  //  act_path_coeff = get_quark_path_coeff();

  ferm_epsilon = weight;

  /* Path coefficients times fermion epsilon */

  c->one_link     = act_path_coeff[0]*ferm_epsilon ;
  c->naik         = act_path_coeff[1]*ferm_epsilon ;
  c->three_staple = act_path_coeff[2]*ferm_epsilon ;
  c->five_staple  = act_path_coeff[3]*ferm_epsilon ;
  c->seven_staple = act_path_coeff[4]*ferm_epsilon ;
  c->lepage       = act_path_coeff[5]*ferm_epsilon ;
}



/*** end of level 3 interface code ****/
/**************************************/


namespace Chroma 
{

  //! Constructor
  /*!
  // Keeping the same interface as for the ordinary staggered 
  // qprop...
  //
  // But the M_ and A_ linop handles are no longer kept
  // (are ignored) -- is there a nice way around this ? 
  // Perhaps not
  */

    // internal variables to this file
static  QDP_ColorVector *out ;
static  QDP_ColorVector *in ; 


  static   QOP_info_t info; // wot for ????
  static   QOP_FermionLinksAsqtad *fla ;  

  // -----------------------------------------------------

  AsqtadCPSWrapperQprop::AsqtadCPSWrapperQprop(
	const	  EvenOddStaggeredTypeFermAct<LatticeStaggeredFermion, 
       	          multi1d<LatticeColorMatrix>, 
                  multi1d<LatticeColorMatrix> >& S_,
	Handle< FermState<LatticeStaggeredFermion,P,Q> > state_,
	  const SysSolverCGParams& invParam_) :
    //    invParam(invParam_),Mass(S_.getQuarkMass()), state(state_)
    invParam(invParam_),Mass(S_.getQuarkMass()), 
state(state_.cast<AsqtadConnectStateBase>())
  {
    // Here is how to get at the gauge links: (Thanks Balint).
    const multi1d<LatticeColorMatrix>& u = state->getLinks();
    const multi1d<LatticeColorMatrix>& u_fat = state->getFatLinks();

#if 0
    multi1d<LatticeColorMatrix> u_with_phases(Nd);
    state.getFermBC().modify(u_with_phases);
#endif

    // add staggered phases to the chroma gauge field
    multi1d<LatticeColorMatrix> u_chroma(Nd);
    // alpha comes from the StagPhases:: namespace
    for(int i = 0; i < Nd; i++) {
      u_chroma[i] = u[i] ; 
      u_chroma[i] *= StagPhases::alpha(i);
    }


  // setup the qdp 
  QOP_asqtad_coeffs_t coeffs;
  load_qop_asqtad_coeffs(&coeffs, 1.0) ;
  Real u0 = 1.0 ; 

  // ---- create the fat links ---------
  QDP_ColorMatrix **uqdp ;
  uqdp = QDP_create_gauge_from_chroma (u_chroma) ;

  QOP_GaugeField *gf;
  gf = QOP_convert_G_from_qdp(uqdp);

  // create the fat links 
  fla = QOP_asqtad_create_L_from_G(&info, &coeffs, gf);

  // remove the space for the gauge configuration
  // gf contains pointers to 4 components of uqdp,
  // so uqdp[i] do not need to be destroyed as well
  QOP_destroy_G(gf) ;
  free(uqdp) ;

}

  
  //! Destructor may not well be automatic
  AsqtadCPSWrapperQprop::~AsqtadCPSWrapperQprop() 
  {
    QOP_asqtad_destroy_L(fla) ;

  }

  SystemSolverResults_t
  AsqtadCPSWrapperQprop::operator() (LatticeStaggeredFermion& psi, 
				     const LatticeStaggeredFermion& q_source) const
  {
    SystemSolverResults_t res;

    QDPIO::cout << "(Hardwired) Mass is " << Mass << endl;

    out = QDP_create_V();
    in = QDP_create_V();


  // convert "q_source" to in (qdp format)
  convert_chroma_to_qdp(in,q_source) ;
  convert_chroma_to_qdp(out,psi) ;

  QOP_ColorVector *qopout, *qopin;
  qopout = QOP_create_V_from_qdp(out);
  qopin = QOP_create_V_from_qdp(in);

  QOP_invert_arg_t inv_arg;
  QOP_resid_arg_t res_arg;
  res_arg.rsqmin = toFloat(invParam.RsdCG)*
    toFloat(invParam.RsdCG);
  inv_arg.max_iter = invParam.MaxCG;
  inv_arg.restart = 1 ;
  inv_arg.evenodd = QOP_EVEN;


  QLA_Real mass = toFloat(Mass) ;
  printf("level3 asqtad inverter mass = %f\n",mass) ; 

  // invert
  QOP_asqtad_invert_all(&info, fla, &inv_arg, &res_arg, mass, qopout, qopin);

  cout << "QOP Inversion results\n" ;
  printf("QOP iters = %d\n",res_arg.final_iter) ;
  printf("QOP ||r||^2 = %g\n",res_arg.final_rsq) ;
  res.n_count = res_arg.final_iter ;
  res.resid = res_arg.final_rsq; // check norm convention

  // convert out to psi (qdp++ format)
  QOP_extract_V_to_qdp(out,qopout) ;
  convert_qdp_to_chroma(psi,out) ;



  QOP_destroy_V(qopout);
  QOP_destroy_V(qopin);

  QDP_destroy_V(out);
  QDP_destroy_V(in);

  return res;
  }

}


void setup_levelthree(multi1d<int> nrow, int *argc, char ***argv )
{
  //  QDP_start(argc, argv);
  int tmp = ::QDP_initialize(argc,argv  ) ;

  printf("Setting up the level 3 code\n");
  printf("----------------------------\n");

  int ndim = 4;
  int *lattice_size;
  lattice_size = (int *) malloc(ndim*sizeof(int));
  for(int i=0 ; i < ndim ; ++i)
    lattice_size[i] = nrow[i] ;

  printf("QDP/QOP lattice volume = %d %d %d %d\n",
	 lattice_size[0],lattice_size[1],
	 lattice_size[2],lattice_size[3] ) ; 

   QDP_set_latsize(ndim, lattice_size);
   QDP_create_layout();

  /** set up the level 3 code **/
  QOP_layout_t qoplayout;
  qoplayout.latdim = ndim;
  qoplayout.latsize = (int *) malloc(ndim*sizeof(int));
  for(int i=0; i<ndim; i++) {
    qoplayout.latsize[i] = lattice_size[i];
  }
  qoplayout.machdim = -1;
  QOP_init(&qoplayout);


  free(lattice_size) ;

}


