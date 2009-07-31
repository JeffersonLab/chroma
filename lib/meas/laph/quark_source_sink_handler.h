#ifndef QUARK_SOURCE_SINK_HANDLER_H
#define QUARK_SOURCE_SINK_HANDLER_H

#include "qdp.h"
#include "chromabase.h"
#include "util/ferm/subset_vectors.h"
#include "meas/inline/io/named_objmap.h"
#include <vector>
#include "gauge_configuration_info.h"
#include "gauge_configuration_handler.h"
#include "xml_help.h"
#include "field_smearing_info.h"
#include "field_smearing_handler.h"
#include "dilution_scheme_info.h"
#include "laph_noise_info.h"
#include "quark_action_info.h"
#include "inverter_info.h"
#include "util/ferm/key_val_db.h"

namespace Chroma {
  namespace LaphEnv {


// ****************************************************************
// *                                                              *
// *  "QuarkSourceSinkHandler" handles access to the smeared gauge  *
// *  field and the eigenvectors of the Laplacian used in the     *
// *  quark field Laplacian Heaviside (Laph) smearing.            *
// *  See declarations below for available member functions.      *
// *                                                              *
// *  Synopsis of important tasks:                                *
// *     -- computes smeared gauge field (locally stored)         *
// *     -- computes Laplacian eigenvectors (NamedObjMap stored)  *
// *                                                              *
// *                                                              *
// *  Internally maintains a QuarkSourceSinkInfo.  Requires         *
// *  connection with an externally maintained                    *
// *  GaugeConfigurationHandler.                                  *
// *                                                              *
// ****************************************************************



class QuarkSourceSinkHandler
{


   struct Key    // used for the file database
    {
       LaphNoiseInfo *noise;  // pointer so can have a default constructor
       int source_time;       // 0..Nt-1 for sinks, Nt for source
       int dilution_index;

       Key();
       Key(const LaphNoiseInfo& in_noise, int in_time, int in_dil_ind);
       Key(const Key& in);
       ~Key();
       Key& operator=(const Key& in);

       bool operator<(const Key& rhs) const;
       void binaryWrite(BinaryWriter& out) const;
       void binaryRead(BinaryReader& in);
       
     //  friend void write(BinaryWriter&, const QuarkSourceSinkHandler::Key&);
     //  friend void read(BinaryReader&, QuarkSourceSinkHandler::Key&);

    };


       // pointers to needed sub-handlers (managed by this handler)
 
   GaugeConfigurationHandler* uPtr;
   FieldSmearingHandler* smearPtr;

       // pointers to internal infos (managed by this handler
       // with new and delete)

   DilutionSchemeInfo *dilPtr;
   QuarkActionInfo *qactionPtr;
   InverterInfo *invertPtr;


       // storage and/or references to internal data

   vector<std::string> fileNames;            // files to handle
   int fileMode;                // 0 = read_only, 1 = must_exist, 
                                // 2 = create_if_not_there
   map<Key,int> fileMap;                    // key -> file index
   map<Key,SerialDBData<LatticeFermion>* > m_storage; // storage of source/sinks


       // Prevent copying ... handler which might contain large
       // amounts of data

   QuarkSourceSinkHandler(const QuarkSourceSinkHandler&);
   QuarkSourceSinkHandler& operator=(const QuarkSourceSinkHandler&);



 public:


   QuarkSourceSinkHandler();

   QuarkSourceSinkHandler(XMLReader& xml_info);

   void setInfo(XMLReader& xml_info);

   void setInfoFromFile(XMLReader& xml_info);

   ~QuarkSourceSinkHandler();

   void clear();

   bool isInfoSet() const;

   string outputInfo() const;



   const GaugeConfigurationInfo& getGaugeConfigurationInfo() const 
      {return uPtr->getInfo();}
   const FieldSmearingInfo& getFieldSmearingInfo() const 
      {return smearPtr->getInfo();}
   const DilutionSchemeInfo& getDilutionSchemeInfo() const 
      {return *dilPtr;}
   const QuarkActionInfo& getQuarkActionInfo() const 
      {return *qactionPtr;}
   const InverterInfo& getInverterInfo() const 
      {return *invertPtr;}

   int getNumberOfDilutionProjectors() const 
    {return dilPtr->getNumberOfProjectors(smearPtr->getInfo());}


        // compute for all dilution indices and dump out to file
        // specified by "file_index";  compute sources for all source
        // times, sinks for just one source time but all sink times
        // (file_index == -1 means use the **last** file in the list)

   void computeSource(XMLReader& xml_in, int file_index = -1);

   void computeSource(const LaphNoiseInfo& noise, int file_index = -1);

   void computeSink(XMLReader& xml_in, int file_index = -1);

   void computeSink(const LaphNoiseInfo& noise, int source_time, int file_index = -1);


   void eraseSink(XMLReader& xml_in);
   
   void eraseSink(const LaphNoiseInfo& noise, int source_time);
   
   void eraseSource(XMLReader& xml_in);
   
   void eraseSource(const LaphNoiseInfo& noise);
   
   

        // read from file (loop through file list to find) for one
        // particular dilution index

   void setSources(const LaphNoiseInfo& noise, int dilution_index);

   void setSink(const LaphNoiseInfo& noise, int source_time, int dilution_index);



   const LatticeFermion& getSources(const LaphNoiseInfo& noise,
                                    int dilution_index);

   const LatticeFermion& getSink(const LaphNoiseInfo& noise, int source_time,
                                 int dilution_index);



        // remove from internal memory

   void removeSourceData(const LaphNoiseInfo& noise, int dilution_index);

   void removeSinkData(const LaphNoiseInfo& noise, int source_time, int dilution_index);

   void clearData();   




 private:


   void create_handlers();
   void destroy_handlers();
   void setFileNames(XMLReader& record_xml);
   void set_info(XMLReader& xml_info);
   void setup_file_map();

   int file_index_helper(int file_index);
   const multi1d<LatticeColorVector>& set_up_laph_eigenvectors();

   friend void write(BinaryWriter& out, const QuarkSourceSinkHandler::Key& key );
   friend void read(BinaryReader& in, QuarkSourceSinkHandler::Key& key);


};



void write(BinaryWriter& out, const QuarkSourceSinkHandler::Key& key );
void read(BinaryReader& in, QuarkSourceSinkHandler::Key& key);


// ***************************************************************
  }
}
#endif  
