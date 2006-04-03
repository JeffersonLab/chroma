// -*- C++ -*-
// $Id: dilute_gauss_src_s.h,v 3.0 2006-04-03 04:59:06 edwards Exp $
/*! \file
 *  \brief Dilute Gaussian sources
 */


//
//
//

#ifndef  DILUTE_GAUSS_SRC_INC
#define  DILUTE_GAUSS_SRC_INC 

namespace Chroma {

  //! Diluted Gauusian-source
  /*! @ingroup sources */
  void gaussian_on_timeslice(LatticeStaggeredFermion& a, 
			     int slice,
			     int mu);
  //! Diluted Gauusian-source
  /*! @ingroup sources */
  void gaussian_on_parity(LatticeStaggeredFermion& a,
			  int parity);
  //! Diluted Gauusian-source
  /*! @ingroup sources */
  void gaussian_color_src(LatticeStaggeredFermion& a,
			  int color_index);
  //! Diluted Gauusian-source
  /*! @ingroup sources */
  void gaussian_color_src_on_slice(LatticeStaggeredFermion& a, 
				   int color_index,
				   int slice, int mu);
  //! Diluted Gauusian-source
  /*! @ingroup sources */
  void gaussian_color_src_on_parity(LatticeStaggeredFermion& a, 
				    int color_index,
				    int parity);
  //! Diluted Gauusian-source
  /*! @ingroup sources */
  void gaussian_parity_src_on_slice(LatticeStaggeredFermion& a,
				    int parity,
				    int slice,
				    int mu);
  //! Diluted Gauusian-source
  /*! @ingroup sources */
  void gaussian_on_mod_timeslice(LatticeStaggeredFermion& a,
				 int slice,
				 int mu,
				 int seperation);
  //! Diluted Gauusian-source
  /*! @ingroup sources */
  void gaussian_on_corner(LatticeStaggeredFermion& a,
			  int corner_index);

  //! Diluted Gauusian-source
  /*! @ingroup sources */
  void gaussian_corner_on_dbl_slice(LatticeStaggeredFermion& a,
				    int corner_index,
				    int slice, int mu);

  //! Diluted Gauusian-source
  /*! @ingroup sources */
  void gaussian_corner_on_mod_dbl_slice(LatticeStaggeredFermion& a,
					int corner_index,
					int slice, int mu, int seperation);
  //! Diluted Gauusian-source
  /*! @ingroup sources */
  void gaussian_color_src_on_mod_slice(LatticeStaggeredFermion& a,
				       int color_index, int slice, int mu,
				       int seperation);

}  // end namespace Chroma

#endif
