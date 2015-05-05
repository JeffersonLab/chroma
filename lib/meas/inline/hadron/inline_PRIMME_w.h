// Use PRIMME from Andreas to find the low mode eigenspace of MdagM. - Arjun

#ifndef __inline_PRIMME_h__
#define __inline_PRIMME_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"
//Add extra file for prop.
#include "io/qprop_io.h"
//For parscalar version, qmp is needed.
#ifdef USE_QMP
#include <qmp.h>
#define fprintf if(QMP_get_node_number()==0)fprintf
#define printf  if(QMP_get_node_number()==0)printf
#endif
#include "actions/ferm/invert/qop_mg/syssolver_qop_mg_params.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlinePrimmeEnv 
  {
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlinehadron */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);

      unsigned long     frequency;

      struct Param_t
      {	

	//Additional params for primme below.
	ChromaProp_t    prop;
	GroupXML_t      fermact;          /*!< fermion action */
	int    nev;          // Number of eigenvalues to be calculated
  	double tol;      // Relatice accuracy of eigenvalues and vectors 
  	int    printLevel;   // Debug

  	double lambda_C;     // cutoff 
  	double lambda_max;   // max eig 
  	int    n_cheb;       // Rank of T_n(x)
	std::string file_name;    // Filename for the eigenvalues/eigenvectors
	int maxBasisSize;       // Determines the maximum basis size passed into primme, if it's 0 nothing is set and defaults are used.
	std::string primme_method;    // String that contains information about what method primme should call.
	SysSolverQOPMGParams invParam; //Params for Multigrid, only visible if turned on.

      };

      struct NamedObject_t
      {
	std::string     gauge_id;      /*!< Gauge field */
      };
      //New structs for primme added below.
      //Holds the input parameters 

      Param_t        param;      /*!< Parameters */
      NamedObject_t  named_obj;  /*!< Named objects */
      std::string    xml_file;   /*!< Alternate XML file pattern */
    };

    //Global sum declared altough may not be used.
    void globalSumDouble(void *sendBuf, void *recvBuf, int *count, void *params);


    //! Inline task for MdagM eigenvectors 
    /*! \ingroup inlinehadron */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the measurement
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    protected:
      //! Do the measurement
      void func(const unsigned long update_no,
		XMLWriter& xml_out); 

    private:
      Params params;
    };

  } // namespace PrimmeEnv

}

#endif
