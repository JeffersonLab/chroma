// 
/*! \file
 * \brief Fermion Flow driver
 *
 */

#include "inline_gluon_cc.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "util/info/unique_id.h"
#include "util/gauge/key_glue_matelem.h"
#include "util/ft/sftmom.h"

#include "util/ferm/key_val_db.h"
#include "util/ferm/key_hadron_2pt_corr.h"

#include "meas/smear/displace.h"


namespace Chroma 
{ 

    //! read input -- gauge fields
    void read(XMLReader& xml, const std::string& path, InlineGluonCCParams::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_in", input.gauge_in);
    }

    //! write output -- gauge fields
    void write(XMLWriter& xml, const std::string& path, const InlineGluonCCParams::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_in", input.gauge_in);

      pop(xml);
    }

    //! read input
    void read(XMLReader& xml, const std::string& path, InlineGluonCCParams::Param_t& input)
    {
      XMLReader inputtop(xml, path);
      read(inputtop, "zmax", input.zmax);
      read(inputtop, "db", input.db);
      read(inputtop, "k5", input.k5);
      read(inputtop,"zdir",input.zdir);

    }

    //! write output
    void write(XMLWriter& xml, const std::string& path, const InlineGluonCCParams::Param_t& input)
    {
      push(xml, path);

      write(xml, "zmax", input.zmax);
      write(xml, "db", input.db);
      write(xml, "k5", input.k5);
      write(xml, "zdir", input.zdir);
      
      pop(xml);
    }




  namespace InlineGluonCCEnv 
  {
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineGluonCC(InlineGluonCCParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }
      
    const std::string name = "GLUON_CC";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  } // namespace InlineGluonCCEnv 


    //--------------------------------------------------------------------------
    // Param stuff
    InlineGluonCCParams::InlineGluonCCParams() { frequency = 0; }

    InlineGluonCCParams::InlineGluonCCParams(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// Parameters for gPDF
	read(paramtop, "Param", param);


	// Read in the source configuration info
	read(paramtop, "NamedObject", named_obj);

	// Possible alternate XML file pattern
	if (paramtop.count("xml_file") != 0) 
	{
	  read(paramtop, "xml_file", xml_file);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
	QDP_abort(1);
      }
    }

  void
  InlineGluonCCParams::writeXML(XMLWriter& xml_out, const std::string& path)
  {
    push(xml_out, path);

    write(xml_out, "Param", param);
    write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }




  // The meat of the following function is copied from the topological charge calculation in chroma/lib/meas/glue/qnaive.cc
  void InlineGluonCC::field_str(multi2d<LatticeColorMatrix>& f, multi1d<LatticeColorMatrix>& u)
  {

    /* Local Variables */
    LatticeColorMatrix u_clov_1;
    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;

    Real k1,k2,k3,k4,kk5;

    Real k5=params.param.k5;

    if( Nd != 4 )
      QDP_error_exit("Nd for the topological charge has to be 4 but: ", Nd);


 // Highly improved FST
/*
    k1 = 19.0/9.0 - 55.0 * k5;
    k1 *= 2.0;
    k2 = 1.0/36.0 - 16.0 * k5;
    k2 *= 2.0;
    k3 = 64.0 * k5 - 32.0/45.0;
    k4 = 1.0/15.0 - 6.0 * k5;
    kk5 = k5;
    kk5 *= 2.0;
*/


 // Simple FST

    k1 = 0.125;
    k2 = 0;
    k3 = 0;
    k4 = 0;
    kk5 = 0;

    for( int mu1 =0; mu1<Nd; mu1++)
    {
	for(int nu1 = 0; nu1<mu1; nu1++)
	{
	    if( toBool(k1 != 0) ) {
	          /* First "plus-plus" 1x1 */
      tmp_1 = u[mu1] * shift(u[nu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(adj(u[mu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 = k1 * tmp_1;

      /* First "plus-minus" 1x1 */
      tmp_1 = u[mu1] * shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[mu1]), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= k1 * tmp_1;

      /* First "minus-minus" 1x1 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k1 * tmp_1;
     /* First "minus-plus" 1x1 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(u[nu1], BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k1 * tmp_1;

      }
     if( toBool(k2!=0) ) {
      /* First "plus-plus" 2x2 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += k2 * tmp_1;

      /* First "plus-minus" 2x2 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= k2 * tmp_1;
      /* First "minus-minus" 2x2 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k2 * tmp_1;
      /* First "minus-plus" 2x2 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1), FORWARD,nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k2 * tmp_1;

      }

      if( toBool(k3!=0) ) {
      /* First "plus-plus" 2x1 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[mu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += k3 * tmp_1;

      /* First "plus-minus" 2x1 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[mu1]), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= k3 * tmp_1;

      /* First "minus-minus" 2x1 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]),BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu1],BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k3 * tmp_1;

      /* First "minus-plus" 2x1 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k3 * tmp_1;


      /* First "plus-plus" 1x2 */
      tmp_1 = u[mu1] * shift(u[nu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], FORWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += k3 * tmp_1;

      /* First "plus-minus" 1x2 */
      tmp_1 = u[mu1] * shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= k3 * tmp_1;
      /* First "minus-minus" 1x2 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k3 * tmp_1;

      /* First "minus-plus" 1x2 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(u[nu1], BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k3 * tmp_1;

      }
      if( toBool(k4!=0) ) {

     /* First "plus-plus" 3x1 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], FORWARD, mu1), FORWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1), FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[mu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += k4 * tmp_1;

      /* First "plus-minus" 3x1 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], FORWARD, mu1), FORWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[mu1]), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= k4 * tmp_1;

      /* First "minus-minus" 3x1 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k4 * tmp_1;
      /* First "minus-plus" 3x1 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k4 * tmp_1;


      /* First "plus-plus" 1x3 */
      tmp_1 = u[mu1] * shift(u[nu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], FORWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[nu1]), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += k4 * tmp_1;

      /* First "plus-minus" 1x3 */
      tmp_1 = u[mu1] * shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= k4 * tmp_1;
     /* First "minus-minus" 1x3 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]),BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k4 * tmp_1;
     /* First "minus-plus" 1x3 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(u[nu1], BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[nu1]), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k4 * tmp_1;

      }


      if( toBool(kk5!=0) ) {

      /* First "plus-plus" 3x3 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], FORWARD, mu1), FORWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1), FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1), FORWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1), FORWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, mu1), FORWARD, nu1),FORWARD, nu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[nu1]), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += kk5 * tmp_1;
      /* First "plus-minus" 3x3 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], FORWARD, mu1), FORWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), FORWARD, mu1),BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), FORWARD,mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[mu1]), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1
), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= kk5 * tmp_1;
      /* First "minus-minus" 3x3 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += kk5 * tmp_1;

      /* First "minus-plus" 3x3 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[nu1]), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= kk5 * tmp_1;

      }

      tmp_1=1;
      tmp_2 = adj(u_clov_1);
      u_clov_1 -= tmp_2;

// make traceless
      tmp_2 = tmp_1 * trace(u_clov_1)/Nc;
      u_clov_1 -= tmp_2;


// JKNOTES: I am almost certain I am missing a normalization here. If I had to guess, I need to divide by a power of 2, but I'll check later
//
//QDPIO::cout<< "JK NOTES:: CHECK NORMALIZATION OF THE FIELD STRENGTH TENSOR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;	
      // Fill the field strength tensor.
	    f[mu1][nu1] = u_clov_1;
	    f[nu1][mu1] = -u_clov_1;
	}// close mu1
    }// close nu1




 }// CLose field_str





    // Function call
    void 
    InlineGluonCC::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {


      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	std::string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "GluonPDF");
	write(xml_out, "update_no", update_no);
	write(xml_out, "xml_file", xml_file);
	pop(xml_out);

	XMLFileWriter xml(xml_file);
	func(update_no, xml);
      }
      else
      {
	func(update_no, xml_out);
      }
    }


    // Real work done here
    // Create the diluted source and apply Lanczos quarature 
    void 
    InlineGluonCC::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      push(xml_out, "GluonPDF");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << InlineGluonCCEnv::name << ": gluon matrix element for parton distributions " << std::endl;

      proginfo(xml_out);    // Print out basic program info

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Test and grab a reference to the gauge field
      // -- we really need only two gauge fields --
      multi1d<LatticeColorMatrix> u ; 
     
      push(xml_out,"GaugeFieldInfo");
      XMLBufferWriter gauge_xml;
      try
	{
	  u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_in);
	  TheNamedObjMap::Instance().get(params.named_obj.gauge_in).getRecordXML(gauge_xml);
	}
      catch( std::bad_cast ) 
	{
	  QDPIO::cerr << InlineGluonCCEnv::name << ": caught dynamic cast error" << std::endl;
	  QDP_abort(1);
	}
      catch (const std::string& e) 
	{
	  QDPIO::cerr << InlineGluonCCEnv::name << ": std::map call failed: " << e << std::endl;
	  QDP_abort(1);
	}
      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      pop(xml_out);

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "GaugeObservables", u);
     
	// I hard code the time direction to be Nd-1. I don't feel bad about this because I hard coded the z direction as well.
	int j_decay = Nd-1;

      // BEGIN JOB 
      
// JKNOTES: Create the key for future use. 
//          Fill the unused values
/*	KeyGlueElementalOperator_t key;
	key.t_slice=0;
	key.left=0;
	key.right=0;

// Open the database
    BinaryStoreDB<  SerialDBKey< KeyGlueElementalOperator_t >, SerialDBData< multi1d<ComplexD> >  > qdp_db;


   if(! qdp_db.fileExists(params.param.db)){
          XMLBufferWriter file_xml;

          push(file_xml, "DBMetaData");
         write(file_xml, "id", std::string("glueElemOp"));
          write(file_xml, "lattSize", QDP::Layout::lattSize());
          proginfo(file_xml);    // Print out basic program info
          write(file_xml, "Params", params.param);
          write(file_xml, "Config_info", gauge_xml);
          pop(file_xml);

          std::string file_str(file_xml.str());
          qdp_db.setMaxUserInfoLen(file_str.size());

          qdp_db.open(params.param.db, O_RDWR | O_CREAT, 0664);
          qdp_db.insertUserdata(file_str);
     }  // close if qdp_db.fileExists
     else
     {
         qdp_db.open(params.param.db, O_RDWR, 0664);
     } // close if else qdp_db.fileExists 

*/

// ADAT Doesn't like reading the glue elemental keys so I am forced to hack a weird key...........
KeyHadron2PtCorr_t key;
key.num_vecs=0;

key.src_name="glue";
key.src_smear="P";
key.src_spin=0;

key.snk_name="glue";
key.snk_smear="P";
key.snk_spin=0;

BinaryStoreDB<  SerialDBKey< KeyHadron2PtCorr_t >, SerialDBData< multi1d<ComplexD> >  > qdp_db;

  if(! qdp_db.fileExists(params.param.db)){
          XMLBufferWriter file_xml;

          push(file_xml, "DBMetaData");
          write(file_xml, "id", std::string("hadron2Pt"));
          write(file_xml, "lattSize", QDP::Layout::lattSize());
	  write(file_xml, "decay_dir", j_decay);
          proginfo(file_xml);    // Print out basic program info
          write(file_xml, "Params", params.param);
          write(file_xml, "Config_info", gauge_xml);
          pop(file_xml);

          std::string file_str(file_xml.str());
          qdp_db.setMaxUserInfoLen(file_str.size());

          qdp_db.open(params.param.db, O_RDWR | O_CREAT, 0664);
          qdp_db.insertUserdata(file_str);
     }  // close if qdp_db.fileExists
     else
     {
         qdp_db.open(params.param.db, O_RDWR, 0664);
     } // close if else qdp_db.fileExists 



	// Calculate the spatial volume
	int space_vol=1;
	for(int dirs=0; dirs <Nd-1; dirs++)
		space_vol*=Layout::lattSize()[dirs];

	// Create the field strength tensor (fst)
	multi2d<LatticeColorMatrix> f(Nd,Nd);
	field_str(f,u);


//open hacky tests
/*
{
   multi1d<int> loc(Nd);
   loc[0]=0;loc[1]=0;loc[2]=0;loc[3]=0;
   ColorMatrix peeker = peekSite(f[1][2],loc);
   int c1=0; int c2 =0;
   Complex res = peekColor(peeker,c1,c2);
   QDPIO::cout << "fieldStr at c1 "<<c1<< " and c2 " <<c2 << " is "<< res<<std::endl;
}
{
   multi1d<int> loc(Nd);
   loc[0]=0;loc[1]=0;loc[2]=0;loc[3]=0;
   LatticeColorMatrix temporary = shift(f[1][3],FORWARD,2);
   ColorMatrix peeker = peekSite(temporary,loc);
   int c1=0; int c2 =0;
   Complex res = peekColor(peeker,c1,c2);
   QDPIO::cout << "fieldStr at c1 "<<c1<< " and c2 " <<c2 << " is "<< res<<std::endl;

}
*/
//close hacky tests



	
	// Set the z direction
	int z_dir = params.param.zdir;

// Store the value with no separation

      // This is hard coded assuming 4 dimensions
      // Make the F^{tx}F^{ty} - F^{ty}F^{tx}
//      LatticeColorMatrix tempF2 = f[Nd-1][0]*f[Nd-1][1] - f[Nd-1][1]*f[Nd-1][0];
//      LatticeComplex F2= trace(tempF2);
      // Make the F\tilde{F}
      // Raza should check the signs on the below and the factor of 1/2
/*    // This is for F\tildeF
      LatticeColorMatrix tempF2 = f[0][1]*f[2][3];
      tempF2 -= f[0][2]*f[1][3];
      tempF2 += f[0][3]*f[1][2];
      tempF2 += f[1][2]*f[0][3];
      tempF2 -= f[1][3]*f[0][2];
      tempF2 += f[2][3]*f[0][1];
      
      LatticeComplex F2fixed = trace(tempF2);
      LatticeComplex F2 = trace(tempF2);

*/


// If changing don't forget to remove } for the for loops at the bottom marked by a comment

for(int m = 0; m<Nd; m++){
for(int n = 0; n<m; n++){
for(int r = 0; r<Nd; r++){
for(int s = 0; s<r; s++){
for(int m2 = 0; m2<Nd; m2++){
for(int n2 = 0; n2<m2; n2++){
for(int r2 = 0; r2<Nd; r2++){
for(int s2 = 0; s2<r2; s2++){

QDPIO::cout << "Op "<< m << n << r << s << " and " << m2 << n2 << r2 << s2 << std::endl;

     LatticeComplex F2fixed = -trace(f[m][n]*f[r][s]);
     LatticeComplex F2 = -trace(f[m2][n2]*f[r2][s2]);


      // Use the sftmom code to create a set for time slices
       			SftMom phases(0, true, j_decay);


			multi1d<ComplexD> corr = sumMulti(F2fixed*F2,phases.getSet());

			multi1d<int> sep(Nd);
			for(int loop = 0; loop<Nd; loop++) {
				sep[loop]=0;
			} // close for loop
			// mom now holds the FST separations
			key.mom=sep;


                        // These variables are now depreciated.
			multi1d<int>lor(4);
			lor[0]=m; lor[1]=n; lor[2]=r; lor[3]=s;
			key.src_lorentz=lor;
			lor[0]=m2; lor[1]=n2; lor[2]=r2; lor[3]=s2;
			key.snk_lorentz=lor;

		// JK NOTES: Save the data
                  SerialDBKey<KeyHadron2PtCorr_t> K;
                  K.key() = key ;
                  SerialDBData< multi1d<ComplexD> >  V ;
			V.data()=corr;
			qdp_db.insert(K,V);

    for(int fb = -1; fb <=1; fb +=2)
    {
	// Create the copy of the fst to shift
	LatticeComplex G2(F2);

	
	// Loop over all lengths of the separation
	for(int z = 1; z<=params.param.zmax ; z++)
	{

		// The QDP Manual warns against using lines such as a = shift(a)
		// temp will be used for these sort of lines
	 	LatticeComplex temp;

		// Displace the FST
		temp=shift(G2,fb, z_dir);

		G2=temp;

       			SftMom phases(0, true, j_decay);


			multi1d<ComplexD> corr = sumMulti(F2fixed*G2,phases.getSet());



			multi1d<int> sep(Nd);
			for(int loop = 0; loop<Nd; loop++) {
				if(loop==z_dir) sep[loop]=fb*z;
				else sep[loop]=0;
			} // close for loop
			// mom now holds the FST separations
			key.mom=sep;


                        // This variable is now depreciated.
			multi1d<int>lor(4);
			lor[0]=m; lor[1]=n; lor[2]=r; lor[3]=s;
			key.src_lorentz=lor;
			lor[0]=m2; lor[1]=n2; lor[2]=r2; lor[3]=s2;
			key.snk_lorentz=lor;

		// JK NOTES: Save the data
                  SerialDBKey<KeyHadron2PtCorr_t> K;
                  K.key() = key ;
                  SerialDBData< multi1d<ComplexD> >  V ;
			V.data()=corr;
			qdp_db.insert(K,V);

	} // close for z

    } // close for fb
      // END JOB

} // close for m
} // close for n
} // close for r
} // close for s
} // close for m2
} // close for n2
} // close for r2
} // close for s2
      
     qdp_db.close(); 
      pop(xml_out);  // GluonPDF
      
      snoop.stop();
      QDPIO::cout << InlineGluonCCEnv::name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << std::endl;
      
      QDPIO::cout << InlineGluonCCEnv::name << ": ran successfully" << std::endl;
      
      END_CODE();
    }
  

} // namespace Chroma
