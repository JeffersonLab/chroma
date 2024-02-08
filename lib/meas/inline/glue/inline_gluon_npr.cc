/*! \file
 *  \brief Inline gluon NPR with CDER
 *
 * Author: Joe Kaprie
 */

#include "inline_gluon_npr.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "util/gauge/taproj.h"


#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

namespace Chroma 
{ 

  namespace InlineGluonNPREnv 
  { 
    //! GluonNPR input
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

      read(paramtop, "r1", param.r1);
      read(paramtop, "r2", param.r2);
      read(paramtop, "rho", param.rho);
      read(paramtop, "p", param.p);
    }

    //! GluonNPR output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& param)
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);
      write(xml, "r1", param.r1);
      write(xml, "r2", param.r2);
      write(xml, "rho", param.rho);
      write(xml, "p", param.p);

      pop(xml);
    }


    //! GluonNPR input
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "obs_id", input.obs_id);
    }

    //! GluonNPR output
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

      const std::string name = "GLUON_NPR";
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

    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "GluonNPR");
      write(xml_out, "update_no", update_no);
      write(xml_out, "param", params.param);
      write(xml_out, "named_obj", params.named_obj);


      if(params.param.rho < 0 || params.param.rho > Nd-1)
        QDP_error_exit("Value of rho is unacceptable");
      if(params.param.p[params.param.rho] != 0 )
        QDP_error_exit("The momentum should be 0 in the gluon's vector direction");

      // Grab the gauge links
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      // Grab the observable
      const LatticeComplex O =
        TheNamedObjMap::Instance().getData< LatticeComplex >(params.named_obj.obs_id); 

      // Get gluon field in rho direction
      LatticeColorMatrix Arho=u[params.param.rho];
//JK DEBUGING:      taproj(Arho);

      // The lattice coordinates for setting the phase
      multi1d<LatticeInteger> my_coord(Nd);
      for (int mu=0; mu < Nd; ++mu)
        my_coord[mu] = Layout::latticeCoordinate(mu);

      // Complex for storing final result
      Complex Lambda(0);
      

      // Loop over the coordinate of the operator in question
      // I am hardcoding 4d and you can't stop me
      multi1d<int> ndim=Layout::lattSize();
      for(int x=0; x < ndim[0]; x++)
      {
       for(int y=0; y < ndim[1]; y++)
       {
        for(int z=0; z < ndim[2]; z++)
        {
         for(int t=0; t < ndim[3]; t++)
         {
           // Grab operator at coordinate coor
           multi1d<int> coor(Nd);
           coor[0] = x; coor[1] = y; coor[2] = z; coor[3] = t;
           Complex Ox = peekSite(O,coor);

           // A set for locations close to coor
           set_r1.make(CloseRangeFunc(params.param.r1, coor));
           set_r2.make(CloseRangeFunc(params.param.r2, coor));

           // Create phase for gluon field FT
           LatticeReal p_dot_x = 0.;
           const Real twopi = 6.283185307179586476925286;
           for(int m = 0; m < Nd; ++m) {
             p_dot_x += LatticeReal(my_coord[m] - coor[m]) * twopi * Real(params.param.p[m]) / ndim[m];
           }
           LatticeComplex phase = cmplx(cos(p_dot_x), sin(p_dot_x)) ;

           // FT the gluon fields for +- p 
           ColorMatrix Ap = sum(Arho * phase, set_r1[0]) ;
           ColorMatrix Amp = sum(Arho * conj(phase), set_r2[0]);

           // Sum up coorelation function
           Lambda += Ox * trace(Ap*Amp) ;
         
         }//close for t
        }//close for z
       }//close for y
      }//close for x

      write(xml_out, "C3", Lambda);

      // Now for the gluon 2pt
      // Create phase for gluon field FT
      LatticeReal p_dot_x = 0.;
      const Real twopi = 6.283185307179586476925286;
      for(int m = 0; m < Nd; ++m) {
        p_dot_x += LatticeReal(my_coord[m]) * twopi * Real(params.param.p[m]) / ndim[m];
      }
      LatticeComplex phase = cmplx(cos(p_dot_x), sin(p_dot_x)) ;

      // FT the gluon fields for +- p
      ColorMatrix Ap = sum(Arho * phase);
      ColorMatrix Amp = sum(Arho * conj(phase));
      Complex C2 = trace(Ap * Amp);

      write(xml_out, "C2", C2);


      pop(xml_out); // pop("GluonNPR");

      END_CODE();
    } 

  }

}
