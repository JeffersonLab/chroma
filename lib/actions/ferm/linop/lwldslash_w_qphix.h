// -*- C++ -*-
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#ifndef __lwldslash_qphix_h__
#define __lwldslash_qphix_h__

#include "qphix_singleton.h"

#ifdef CHROMA_BUILDING_QPHIX_DSLASH
#warning USING QPHIX DSLASH
#undef DEBUG_QPHIX_DSLASH


#include "chroma_config.h"
#include "actions/ferm/linop/lwldslash_base_w.h"
#include "io/aniso_io.h"
#include "state.h"
#include "qdp_datalayout.h"
#include "qphix/geometry.h"
#include "qphix/dslash_def.h"
#include "qphix/qdp_packer.h"
#include "actions/ferm/invert/qphix/qphix_vec_traits.h"



#ifdef S
warning "S is defined"
#endif

namespace Chroma 
{ 
  //! General Wilson-Dirac dslash
  /*!
   * \ingroup linop
   *
   * DSLASH
   *
   * This routine is specific to Wilson fermions!
   *
   * Description:
   *
   * This routine applies the operator D' to Psi, putting the result in Chi.
   *
   *	       Nd-1
   *	       ---
   *	       \
   *   chi(x)  :=  >  U  (x) (1 - isign gamma  ) psi(x+mu)
   *	       /    mu			  mu
   *	       ---
   *	       mu=0
   *
   *	             Nd-1
   *	             ---
   *	             \    +
   *                +    >  U  (x-mu) (1 + isign gamma  ) psi(x-mu)
   *	             /    mu			   mu
   *	             ---
   *	             mu=0
   *
   */
  template<typename T,typename P, typename Q>
  class QPhiXWilsonDslash : public WilsonDslashBase<T,P,Q>
  {
  public:
	  // Floating Type <float or double>
	    using FT = typename WordType<T>::Type_t;

	    // Save fingers...
	    template<typename TT>
	   		using VecTraits = Chroma::QPhiXVecTraits::VecTraits<TT>;

	    // Geometry Type
	    using Geom = typename QPhiX::Geometry<FT,VecTraits<FT>::Vec, VecTraits<FT>::Soa, false>;

	    // Dslash Type
	    // FIXME: IT is possible that some noddy user applies this to a non SU(3) field. It should still
	    //   Give the same result as a regular Dslash, so I will for safety set the compress to false here.
	    using Dsl = typename QPhiX::Dslash<FT,VecTraits<FT>::Vec, VecTraits<FT>::Soa, false>;

	    // Spinor and Gauge Types
	    using Spinor = typename Geom::FourSpinorBlock;
	    using Gauge = typename Geom::SU3MatrixBlock;



  public:
    //! Empty constructor. Must use create later

	    // Originally had constructor args:     : theGeom(nullptr), theDslash(nullptr), packed_gauge{nullptr,nullptr}
	    	// But now these are static
	QPhiXWilsonDslash()
#ifndef CHROMA_STATIC_PACKED_GAUGE
  	  : packed_gauge{nullptr,nullptr}
#endif
   {}




    //! Full constructor
	// Originally had constructor args:     : theGeom(nullptr), theDslash(nullptr), packed_gauge{nullptr,nullptr}
		// But now these are static
    QPhiXWilsonDslash(Handle< FermState<T,P,Q> > state)
#ifndef CHROMA_STATIC_PACKED_GAUGE
  	  : packed_gauge{nullptr,nullptr}
#endif
    {
    	create(state);
    }

    //! Full constructor with anisotropy
    // Originally had constructor args:     : theGeom(nullptr), theDslash(nullptr), packed_gauge{nullptr,nullptr}
    	// But now these are static
    QPhiXWilsonDslash(Handle< FermState<T,P,Q> > state,
		    const AnisoParam_t& aniso_)
#ifndef CHROMA_STATIC_PACKED_GAUGE
  	  : packed_gauge{nullptr,nullptr}
#endif
    {
    	create(state,aniso_);
    }

    //! Full constructor with general coefficients
    // Originally had constructor args:     : theGeom(nullptr), theDslash(nullptr), packed_gauge{nullptr,nullptr}
    	// But now these are static
    QPhiXWilsonDslash(Handle< FermState<T,P,Q> > state,
		    const multi1d<Real>& coeffs_)
#ifndef CHROMA_STATIC_PACKED_GAUGE
  	  : packed_gauge{nullptr,nullptr}
#endif
    {
    	create(state,coeffs_);
    }

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > state)
    {
      multi1d<Real> cf(Nd);
      cf = 1.0;
      create(state, cf);
    }
    //! Creation routine with anisotropy
    void create(Handle< FermState<T,P,Q> > state, 
		const AnisoParam_t& aniso_) {
    	  START_CODE();

    	    create(state, makeFermCoeffs(aniso_));
    	    END_CODE();
    }

    //! Full constructor with general coefficients
    void create(Handle< FermState<T,P,Q> > state, const multi1d<Real>& coeffs_) {
    	START_CODE();


    	// Save a copy of the fermbc
    	fbc = state->getFermBC();
    	coeffs = coeffs_;

    	// Sanity check
    	if (fbc.operator->() == 0) {
    		QDPIO::cerr << "QPhiXWilsonDslash: error: fbc is null" << std::endl;
    		QDP_abort(1);
    	}

    	// Get the geometry:
    	QPhiX::QPhiXCLIArgs& QPhiXParams  = TheQPhiXParams::Instance();



    	if( theGeom == nullptr ) {
#ifdef DEBUG_QPHIX_DSLASH
    		QDPIO::cout << "Allocating Geometry" <<std::endl;
    		QDPIO::cout << "   By="<< QPhiXParams.getBy() << std::endl;
    		QDPIO::cout << "   Bz="<< QPhiXParams.getBz() << std::endl;
    		QDPIO::cout << "   NCores="<< QPhiXParams.getNCores() << std::endl;
    		QDPIO::cout << "   Sy="<< QPhiXParams.getSy() << std::endl;
    		QDPIO::cout << "   Sz="<< QPhiXParams.getSz() << std::endl;
    		QDPIO::cout << "   MinCt="<< QPhiXParams.getMinCt() << std::endl;

#endif
    		theGeom = new Geom(Layout::subgridLattSize().slice(),
    						   QPhiXParams.getBy(),
							   QPhiXParams.getBz(),
							   QPhiXParams.getNCores(),
							   QPhiXParams.getSy(),
							   QPhiXParams.getSz(),
							   0,0,QPhiXParams.getMinCt()); // Set Pads to zero explicitly

    		if (theGeom == nullptr) {
    			QDPIO::cerr <<"Failed to allocate geometry. Aborting" <<std::endl;
    			QDP_abort(1);
    		}
    	}

    	if( theDslash == nullptr ) {
    		double aniso_fac_s = toDouble(coeffs_[0]);
    		double aniso_fac_t = toDouble(coeffs_[Nd-1]);
    		double t_boundary = toDouble(1); // FBC will apply boundaries...
#ifdef DEBUG_QPHIX_DSLASH
    		QDPIO::cout << "Allocating Dslash: fac_s="<<aniso_fac_s<<" fac_t="<<aniso_fac_t << " t_boundary="<<t_boundary <<std::endl;
#endif
    		theDslash = new Dsl(theGeom,t_boundary,aniso_fac_s, aniso_fac_t);
    		if( theDslash == nullptr ) {
    			QDPIO::cerr << "Failed to allocate Dslash... Aborting " << std::endl;
    			QDP_abort(1);
    		}
    	}

    	if( packed_gauge[0] == nullptr) {
#ifdef DEBUG_QPHIX_DSLASH
    		QDPIO::cout << "Allocating packed gauge field (cb=0)" << std::endl;
#endif

    		packed_gauge[0] = (Gauge*)theGeom->allocCBGauge();
    		if( packed_gauge[0] == nullptr ) {
    			QDPIO::cerr << "Failed to allocate packed gauge_field(cb=0)" << std::endl;
    			QDP_abort(1);
    		}
    	}

    	if( packed_gauge[1] == nullptr) {
#ifdef DEBUG_QPHIX_DSLASH
    	QDPIO::cout << "Allocating packed gauge field (cb=1)" << std::endl;
#endif
    		packed_gauge[1] = (Gauge*)theGeom->allocCBGauge();
    		if( packed_gauge[0] == nullptr ) {
    			QDPIO::cerr << "Failed to allocate packed gauge_field(cb=1)" << std::endl;
    			QDP_abort(1);
    		}
    	}

    	const Q& u = state->getLinks();
#ifdef DEBUG_QPHIX_DSLASH
    	QDPIO::cout << "Packing Gauge Field" << std::endl;
#endif
    	QPhiX::qdp_pack_gauge<>(u,packed_gauge[0],packed_gauge[1],(*theGeom));
    	END_CODE();
    }

    //! No real need for cleanup here
    ~QPhiXWilsonDslash() {

#ifndef CHROMA_STATIC_PACKED_GAUGE
    	// Keeping the packed gauges static shares them betwen objects which may be
    	// dangerous. SO this needs to be explicitly enabled, otherwise we free and reallocate
    theGeom->free(packed_gauge[0]);  packed_gauge[0] = nullptr;
   	theGeom->free(packed_gauge[1]);  packed_gauge[1] = nullptr;
#endif

   	// Never free these.
  //  	delete theDslash; theDslash = nullptr;
  //	delete theGeom; theGeom = nullptr;
    }

    /**
     * Apply a dslash
     *
     * \param chi     result                                      (Write)
     * \param psi     source                                      (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of OUTPUT std::vector               (Read) 
     *
     * \return The output of applying dslash on psi
     */
    void apply(T& chi, const T& psi, 
	       enum PlusMinus isign, int cb) const
    {
    	int source_cb = 1 - cb;
    	int target_cb = cb;
    	int qphix_isign = (isign == PLUS) ? +1 : -1;
    	int cbsize_in_blocks = rb[0].numSiteTable()/VecTraits<FT>::Soa;

    	// Layout etc has been checked on construction and there is no padding
    	Spinor* chi_targ = (Spinor *)(chi.getFjit())+target_cb*cbsize_in_blocks;
    	Spinor* psi_src  = (Spinor *)(psi.getFjit())+source_cb*cbsize_in_blocks;
    	theDslash->dslash(chi_targ,
    	      		 	 psi_src,
						 packed_gauge[target_cb],
						 qphix_isign,
						 target_cb);

    	// Apply Boundary
    	getFermBC().modifyF(chi, QDP::rb[cb]);
    }

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

  protected:
    //! Get the anisotropy parameters
    const multi1d<Real>& getCoeffs() const {return coeffs;}

  private:
    multi1d<Real> coeffs;  /*!< Nd array of coefficients of terms in the action */
    Handle< FermBC<T,P,Q> > fbc;


    // Pointers to held geom and dslash
    static Geom *theGeom;
    static Dsl *theDslash;

    // Pointers to packed gauge field
#ifdef CHROMA_STATIC_PACKED_GAUGE
    static Gauge* packed_gauge[2];
#else
    Gauge* packed_gauge[2];
#endif
  };


  // The Concrete Instances
  using QPhiXWilsonDslashFloating = QPhiXWilsonDslash<LatticeFermion,multi1d<LatticeColorMatrix>,multi1d<LatticeColorMatrix>>;
  using QPhiXWilsonDslashF = QPhiXWilsonDslash<LatticeFermionF,multi1d<LatticeColorMatrixF>,multi1d<LatticeColorMatrixF>>;
  using QPhiXWilsonDslashD = QPhiXWilsonDslash<LatticeFermionD,multi1d<LatticeColorMatrixD>,multi1d<LatticeColorMatrixD>>;


  QPhiXWilsonDslashF::Geom* QPhiXWilsonDslashF::theGeom = nullptr;
  QPhiXWilsonDslashF::Dsl* QPhiXWilsonDslashF::theDslash = nullptr;

#ifdef CHROMA_STATIC_PACKED_GAUGE
  QPhiXWilsonDslashF::Gauge* QPhiXWilsonDslashF::packed_gauge[2] = { nullptr,nullptr };
#endif

  QPhiXWilsonDslashD::Geom* QPhiXWilsonDslashD::theGeom = nullptr;
  QPhiXWilsonDslashD::Dsl* QPhiXWilsonDslashD::theDslash = nullptr;

#ifdef CHROMA_STATIC_PACKED_GAUGE
  QPhiXWilsonDslashD::Gauge* QPhiXWilsonDslashD::packed_gauge[2] = {nullptr,nullptr};
#endif

} // End Namespace Chroma

#endif // Ifdef BUILDING_CHROMA_QPHIX_DSLASH
#endif
