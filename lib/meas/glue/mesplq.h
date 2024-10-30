// -*- C++ -*-

#ifndef __mesplq_h__
#define __mesplq_h__

namespace Chroma 
{

  //! Return the value of the average plaquette normalized to 1
  /*!
   * \ingroup glue
   *
   * \param u         gauge field (Read)
   * \param w_plaq    plaquette average (Write)
   * \param s_plaq    space-like plaquette average (Write)
   * \param t_plaq    time-like plaquette average (Write)
   * \param link      space-time average link (Write)
   */
#if ! defined (QDP_IS_QDPJIT2)
  void MesPlq(const multi1d<LatticeColorMatrixF3>& u, 
	      Double& w_plaq, Double& s_plaq, Double& t_plaq, Double& link);

  void MesPlq(const multi1d<LatticeColorMatrixD3>& u, 
	      Double& w_plaq, Double& s_plaq, Double& t_plaq, Double& link);
#else
  void MesPlq(const multi1d<LatticeColorMatrix>& u, 
	      Double& w_plaq, Double& s_plaq, Double& t_plaq, Double& link);
#endif
  
  //! Return the value of the average plaquette normalized to 1
  /*!
   * \ingroup glue
   *
   * \param u           gauge field (Read)
   * \param w_plaq      plaquette average (Write)
   * \param s_plaq      space-like plaquette average (Write)
   * \param t_plaq      time-like plaquette average (Write)
   * \param plane_plaq  plane plaquette average (Write)
   * \param link        space-time average link (Write)
   */

#if ! defined (QDP_IS_QDPJIT2)
  void MesPlq(const multi1d<LatticeColorMatrixF3>& u, 
	      Double& w_plaq, Double& s_plaq, Double& t_plaq, 
	      multi2d<Double>& plane_plaq,
	      Double& link);

  void MesPlq(const multi1d<LatticeColorMatrixD3>& u, 
	      Double& w_plaq, Double& s_plaq, Double& t_plaq, 
	      multi2d<Double>& plane_plaq,
	      Double& link);
#else
  void MesPlq(const multi1d<LatticeColorMatrix>& u, 
	      Double& w_plaq, Double& s_plaq, Double& t_plaq, 
	      multi2d<Double>& plane_plaq,
	      Double& link);
#endif
  
  //! Print the value of the average plaquette normalized to 1
  /*!
   * \ingroup glue
   *
   * \param xml    plaquette average (Write)
   * \param u      gauge field (Read)
   */
#if ! defined (QDP_IS_QDPJIT2)
  void MesPlq(XMLWriter& xml,
	      const std::string& xml_group,
	      const multi1d<LatticeColorMatrixF3>& u);

  void MesPlq(XMLWriter& xml,
	      const std::string& xml_group,
	      const multi1d<LatticeColorMatrixD3>& u);
#else
  void MesPlq(XMLWriter& xml,
	      const std::string& xml_group,
	      const multi1d<LatticeColorMatrix>& u);
#endif
  
}  // end namespace Chroma

#endif
