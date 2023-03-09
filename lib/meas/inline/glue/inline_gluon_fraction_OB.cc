/*! \file
 *  \brief Inline gluon momentum fraction operator OB
 *  O_mn = 2 Tr[F_ms F_ns]
 *  O_B = O_44 - 1/3 O_ii
 *
 * Author: Joe Karpie
 */

#include "inline_gluon_fraction_OB.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"

#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"


namespace Chroma 
{ 

  namespace InlineGluonMomFracOBEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path)
      {
	Params p(xml_in, path);
	return new InlineMeas(p);
      }

      //! Local registration flag
      bool registered = false;

      const std::string name = "GLUON_MOM_FRAC_OB";
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= CreateGaugeStateEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  

    //! GluonMomFracOB input
    void read(XMLReader& xml, const std::string& path, Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	break;

      default:
	QDPIO::cerr << "Params::Param_t: " << version 
		    << " unsupported." << std::endl;
	QDP_abort(1);
      }
    }

    //! GluonMomFracOB output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& param)
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);

      pop(xml);
    }


    //! GluonMomFracOB input
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "obs_id", input.obs_id);
    }

    //! GluonMomFracOB output
    void write(XMLWriter& xml, const std::string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "obs_id", input.obs_id);

      pop(xml);
    }


    // Params
    Params::Params()
    { 
      frequency = 0; 
    }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// Params
	read(paramtop, "Param", param);

	// Ids
	read(paramtop, "NamedObject", named_obj);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << "Caught Exception reading XML: " << e << std::endl;
	QDP_abort(1);
      }
    }

  // The meat of the following function is copied from the topological charge calculation in chroma/lib/meas/glue/qnaive.cc
  void field_str(multi2d<LatticeColorMatrix>& f, const multi1d<LatticeColorMatrix>& u)
  {

    /* Local Variables */
    LatticeColorMatrix u_clov_1;
    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;

    Real k1,k2,k3,k4,kk5;


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

    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      // Grab the object
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      push(xml_out, "GluonMomFracOB");
      write(xml_out, "update_no", update_no);

      // Create the field strength tensor (fst)
      multi2d<LatticeColorMatrix> F(Nd,Nd);
      field_str(F,u);

      //  O_mn = 2 Tr[F_ms F_ns]
      //  O_B = O_44 - 1/3 O_ii

      LatticeComplex Ob = Real(4./3.) * trace(F[3][0] * F[3][0] + F[3][1] * F[3][1] + F[3][2] * F[3][2]);

      Ob -= Real(4./3.) * trace(F[0][1]*F[0][1] + F[0][2]*F[0][2] + F[1][2]*F[1][2]);

       // Now store the operator to a memory object
      {
        XMLBufferWriter file_xml, record_xml;
        push(file_xml, "Gluon_mom_frac_Ob");
        write(file_xml, "id", int(0));
        write(file_xml, "GluonMomFracOBParams", params.param);
        pop(file_xml);

        // Store the gauge field
        TheNamedObjMap::Instance().create< LatticeComplex >(params.named_obj.obs_id);
        TheNamedObjMap::Instance().getData< LatticeComplex >(params.named_obj.obs_id) = Ob;
        TheNamedObjMap::Instance().get(params.named_obj.obs_id).setFileXML(file_xml);
        TheNamedObjMap::Instance().get(params.named_obj.obs_id).setRecordXML(record_xml);
      }

      pop(xml_out); // pop("GluonMomFracOB");

      END_CODE();
    } 

  }

}
