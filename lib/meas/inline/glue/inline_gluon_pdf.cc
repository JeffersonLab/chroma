// 
/*! \file
 * \brief Fermion Flow driver
 *
 */

#include "inline_gluon_pdf.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "util/info/unique_id.h"
#include "util/gauge/key_glue_matelem.h"
#include "util/ft/sftmom.h"

#include "util/ferm/key_val_db.h"
#include "util/ferm/unsmeared_fsq_disco.h"

#include "meas/smear/displace.h"


namespace Chroma 
{ 

    //! read input -- gauge fields
    void read(XMLReader& xml, const std::string& path, InlineGluonPDFParams::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_in", input.gauge_in);
    }

    //! write output -- gauge fields
    void write(XMLWriter& xml, const std::string& path, const InlineGluonPDFParams::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_in", input.gauge_in);

      pop(xml);
    }

    //! read input
    void read(XMLReader& xml, const std::string& path, InlineGluonPDFParams::Param_t& input)
    {
      XMLReader inputtop(xml, path);
      read(inputtop, "zmax", input.zmax);
      read(inputtop, "db", input.db);
      read(inputtop, "k5", input.k5);
      read(inputtop,"zdir",input.zdir);
      read(inputtop,"smear",input.smear);
      read(inputtop,"mom",input.mom);

    }

    //! write output
    void write(XMLWriter& xml, const std::string& path, const InlineGluonPDFParams::Param_t& input)
    {
      push(xml, path);

      write(xml, "zmax", input.zmax);
      write(xml, "db", input.db);
      write(xml, "k5", input.k5);
      write(xml, "zdir", input.zdir);
      write(xml, "mom", input.mom);
      
      pop(xml);
    }




  namespace InlineGluonPDFEnv 
  {
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineGluonPDF(InlineGluonPDFParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }
      
    const std::string name = "GLUON_PDF";

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
  } // namespace InlineGluonPDFEnv 


    //--------------------------------------------------------------------------
    // Param stuff
    InlineGluonPDFParams::InlineGluonPDFParams() { frequency = 0; }

    InlineGluonPDFParams::InlineGluonPDFParams(XMLReader& xml_in, const std::string& path) 
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
  InlineGluonPDFParams::writeXML(XMLWriter& xml_out, const std::string& path)
  {
    push(xml_out, path);

    write(xml_out, "Param", param);
    write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }




  // The meat of the following function is copied from the topological charge calculation in chroma/lib/meas/glue/qnaive.cc
  void InlineGluonPDF::field_str(multi2d<LatticeColorMatrix>& f, multi1d<LatticeColorMatrix>& u)
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


      // Fill the field strength tensor.
	    f[mu1][nu1] = u_clov_1;
	    f[nu1][mu1] = -u_clov_1;
	}// close mu1
    }// close nu1




 }// CLose field_str





    // Function call
    void 
    InlineGluonPDF::operator()(unsigned long update_no,
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
    InlineGluonPDF::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      push(xml_out, "GluonPDF");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << InlineGluonPDFEnv::name << ": gluon matrix element for parton distributions " << std::endl;

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
	  QDPIO::cerr << InlineGluonPDFEnv::name << ": caught dynamic cast error" << std::endl;
	  QDP_abort(1);
	}
      catch (const std::string& e) 
	{
	  QDPIO::cerr << InlineGluonPDFEnv::name << ": std::map call failed: " << e << std::endl;
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

// ADAT Doesn't like reading the glue elemental keys so I am forced to hack a weird key...........
KeyFSqDiscoOperator_t key;
key.smear = params.param.smear;


//Generating phase for mom transfer
multi1d<int> mom = params.param.mom;
    LatticeComplex qq=1.0 ;
    //IF YOU DO NOT USE A MOMENTUM FACTOR in time just
    // give 3 indices for the momentum 
    for(int k(0);k<mom.size();k++){
      LatticeReal pk = (Layout::latticeCoordinate(k) )
        * mom[k] * twopi / Real(Layout::lattSize()[k]);
      qq = qq*cmplx(cos(pk),sin(pk));
    }

key.mom=mom;


BinaryStoreDB<  SerialDBKey< KeyFSqDiscoOperator_t >, SerialDBData< multi1d<ComplexD> >  > qdp_db;

  if(! qdp_db.fileExists(params.param.db)){
          XMLBufferWriter file_xml;

          push(file_xml, "DBMetaData");
          write(file_xml, "id", std::string("FsqDiscoBlocks"));
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

	// Set the z direction
	int z_dir = params.param.zdir;

// Store the value with no separation

	for( int m =1; m<Nd; m++)
	{
	    for(int n = 0 ; n<m; n++ )
	    {
		for(int r = 1; r<Nd; r++){
		for(int s = 0; s<r; s++){

      // Use the sftmom code to create a set for time slices
       			SftMom phases(0, true, j_decay);


			multi1d<ComplexD> corr = sumMulti(-qq*trace(f[m][n]*f[r][s]),phases.getSet());

			std::vector<int> disp_list;
			key.disp_list = disp_list;


			std::vector<int> left_lorentz = {m,n};
			key.left_lorentz=left_lorentz;
			std::vector<int> right_lorentz = {r,s};
			key.right_lorentz=right_lorentz;

		// JK NOTES: Save the data
                  SerialDBKey<KeyFSqDiscoOperator_t> K;
                  K.key() = key ;
                  SerialDBData< multi1d<ComplexD> >  V ;
			V.data()=corr;
			qdp_db.insert(K,V);

		} // close for s
		} // close for r
	    } //close for n
	}// close for m

    for(int fb = -1; fb <=1; fb +=2)
    {
	// Create the copy of the fst to shift
	multi2d<LatticeColorMatrix> g(f);
//	multi2d<LatticeColorMatrix> h(f);

	// Create the copy of the u_z links to reconnect the FSTs 
	// Needed for adjoint rep definition of gluon PDF
	LatticeColorMatrix uz=1;

	// Create a copy of the links to follow the shifts performed so that consecutive displace calls will work with correct spatial coordinates
	
        std::vector<int> disp_list;

	// Loop over all lengths of the wilson line
	for(int z = 1; z<=params.param.zmax ; z++)
	{
		// The QDP Manual warns against using lines such as a = shift(a)
		// temp will be used for these sort of lines
	 	LatticeColorMatrix temp;

		// Build the link which will go from z to 0
	        temp = displace(u,uz,fb,z_dir);
	        uz = temp;

		// Displace the FST
		for(int mu=0; mu<Nd; mu++)
		{
			for(int nu=0; nu<Nd; nu++)
			{
				temp = displace(u,g[mu][nu],fb,z_dir);
				g[mu][nu] = temp;
//				temp = shift(h[mu][nu],fb,z_dir);
//				h[mu][nu]=temp;
			} // close for nu
		} // close for mu



        // Update the disp list before saving
        disp_list.push_back(fb*(z_dir));
        key.disp_list = disp_list;



	for(int m =1 ; m<Nd; m++)
	{
	    for(int n = 0 ; n<m; n++ )
	    {
		for(int r = 0; r<Nd; r++){
		for(int s = 0; s<Nd; s++){
		    if(r!=s)
		    {

      // Use the sftmom code to create a set for time slices
       			SftMom phases(0, true, j_decay);


			multi1d<ComplexD> corr =  sumMulti(-qq*trace(f[m][n]*g[r][s]*adj(uz)),phases.getSet());



                        std::vector<int> left_lorentz = {m,n};
                        key.left_lorentz=left_lorentz;
                        std::vector<int> right_lorentz = {r,s};
                        key.right_lorentz=right_lorentz;

		// JK NOTES: Save the data
                  SerialDBKey<KeyFSqDiscoOperator_t> K;
                  K.key() = key ;
                  SerialDBData< multi1d<ComplexD> >  V ;
			V.data()=corr;
			qdp_db.insert(K,V);

		    } // close if r!=s
		} // close for s
		} // close for r
	    } // close for n
	 } // close for m
	} // close for z

    } // close for fb
      // END JOB

      
     qdp_db.close(); 
      pop(xml_out);  // GluonPDF
      
      snoop.stop();
      QDPIO::cout << InlineGluonPDFEnv::name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << std::endl;
      
      QDPIO::cout << InlineGluonPDFEnv::name << ": ran successfully" << std::endl;
      
      END_CODE();
    }
  

} // namespace Chroma
