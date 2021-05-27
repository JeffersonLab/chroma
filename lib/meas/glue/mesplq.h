// -*- C++ -*-

#ifndef __mesplq_h__
#define __mesplq_h__

namespace Chroma 
{

  template<typename Q>
  void MesPlq(const multi1d<Q>& u, 
	      multi2d<Double>& plane_plaq, Double& link);

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

  template<typename Q>
  void MesPlq(const multi1d<Q>& u, 
	      Double& w_plaq, Double& s_plaq, Double& t_plaq, 
	      multi2d<Double>& plane_plaq,
	      Double& link);
  template<typename Q>
  void MesPlq(const multi1d<Q>& u, 
	      Double& w_plaq, Double& s_plaq, Double& t_plaq, Double& link);

  //! Print the value of the average plaquette normalized to 1
  /*!
   * \ingroup glue
   *
   * \param xml    plaquette average (Write)
   * \param u      gauge field (Read)
   */
  template<typename Q>
  void MesPlq(XMLWriter& xml, 
	      const std::string& xml_group,
	      const multi1d<Q>& u);

}  // end namespace Chroma

#endif
