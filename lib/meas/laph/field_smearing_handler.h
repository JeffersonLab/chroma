#ifndef FIELD_SMEARING_HANDLER_H
#define FIELD_SMEARING_HANDLER_H

#include "qdp.h"
#include "chromabase.h"
#include "util/ferm/subset_vectors.h"
#include "meas/inline/io/named_objmap.h"
#include "util/gauge/stout_utils.h"
#include "gauge_configuration_handler.h"
#include "field_smearing_info.h"

namespace Chroma {
  namespace LaphEnv {


// ****************************************************************
// *                                                              *
// *  "FieldSmearingHandler" handles access to the smeared gauge  *
// *  field and the eigenvectors of the Laplacian used in the     *
// *  quark field Laplacian Heaviside (Laph) smearing.            *
// *  See declarations below for available member functions.      *
// *                                                              *
// *  Synopsis of important tasks:                                *
// *     -- computes smeared gauge field (locally stored)         *
// *     -- computes Laplacian eigenvectors (NamedObjMap stored)  *
// *                                                              *
// *                                                              *
// *  Internally maintains a FieldSmearingInfo and a              *
// *  GaugeConfigurationHandler.                                  *
// *                                                              *
// *  All Laph Handlers follow the member naming convention:      *
// *                                                              *
// *    compute....()  to do original computation                 *
// *    set...()       to internally set from file or NamedObjMap *
// *                                                              *
// *    get...()       provides access to results                 *
// *                                                              *
// ****************************************************************



class FieldSmearingHandler
{

       // pointers to other handlers (managed by this handler
       // with new and delete)
 
   GaugeConfigurationHandler* uPtr;

       // pointers to internal infos (managed by this handler
       // with new and delete)

   const FieldSmearingInfo* smearPtr;


       // storage and/or references to internal data

   multi1d<LatticeColorMatrix> usmear;   // storage for smeared gauge field
 
   SubsetVectors<LatticeColorVector>* laph_eigvecs;  // stored in NamedObjMap

   string laph_eigvecs_id;

       // prevent copying ... handler might contain large
       // amounts of data

   FieldSmearingHandler(const FieldSmearingHandler&);
   FieldSmearingHandler& operator=(const FieldSmearingHandler&);



 public:


   FieldSmearingHandler();

   FieldSmearingHandler(XMLReader& xml_in);

   void setInfo(XMLReader& xml_in);

   void setInfo(const std::string& header);

   ~FieldSmearingHandler();

   void clear();  // clears everything in the handler


           // access to the info

   bool isInfoSet() const;

   const FieldSmearingInfo& getFieldSmearingInfo() const;

   const GaugeConfigurationInfo& getGaugeConfigurationInfo() const;

   string outputInfo() const;



           // compute the Laph Eigenvectors and put into NamedObjMap
           // "laph_eigvecs" will point to result

   void computeLaphEigenvectors();

           // set Laph eigenvectors from the NamedObjMap

   void setLaphEigenvectors();

           // provides access to the quark smearing eigenvectors

   const multi1d<Real>& getLaphEigenvalues(int t) const;

   const multi1d<LatticeColorVector>& getLaphEigenvectors() const;

           // clear pointer to Laph eigenvectors, and if "namedobj_erase"
           // is true, also remove from the NamedObjMap

   void clearLaphEigenvectors(bool namedobj_erase = false);

   bool isLaphEigenvectorsSet() const {return (laph_eigvecs!=0);}



           // compute the smeared gauge field and store internally 
           // in "usmear"

   void computeSmearedGaugeField();

           // provides access to the smeared gauge field

   const multi1d<LatticeColorMatrix>& getSmearedGaugeField() const;

           // clear the smeared gauge field from "usmear"

   void clearGaugeField();

   bool isSmearedGaugeFieldSet() const {return (usmear.size()!=0);}





 private:

   void set_info(XMLReader& xml_in);
   void clear_data();
   void create_handlers();
   void destroy_handlers();
   void destroy();
   void check_match(XMLReader& record_xml);


};


// ***************************************************************
  }
}
#endif  
