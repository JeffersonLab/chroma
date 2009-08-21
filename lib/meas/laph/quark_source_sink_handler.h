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
#include "quark_info.h"
#include "inverter_info.h"
#include "util/ferm/key_val_db.h"

namespace Chroma {
  namespace LaphEnv {


// ****************************************************************
// *                                                              *
// *  "QuarkSourceSinkHandler" handles computation of and         *
// *  subsequent access to the quark source and sink functions.   *
// *                                                              *
// *  Required XML input for setting the handler info:            *
// *                                                              *
// *   <QuarkSourceSinkInfo>                                      *
// *      <FileNameList> ...  </FileNameList>  (required)         *
// *      <InvertParam>  ...  </InvertParam>   (required)         *
// *      <StoutLaphSmearing>      ...  </StoutLaphSmearing>      *
// *      <GaugeConfigurationInfo> ...  </GaugeConfigurationInfo> *
// *      <QuarkInfo>              ...  </QuarkInfo>              *
// *      <LaphDilutionScheme>     ...  </LaphDilutionScheme>     *
// *   </QuarkSourceSinkInfo>                                     *
// *                                                              *
// *  The file list and the inverter parameters are mandatory.    *
// *  The remaining info items, if absent, are extracted from     *
// *  the files in file list.  If no files yet exist, then these  *
// *  info items are mandatory.                                   *
// *  All Laph Handlers follow the member naming convention:      *
// *                                                              *
// *    compute....()  to do original computation                 *
// *    set...()       to internally set from file or NamedObjMap *
// *                                                              *
// *    get...()       provides access to results                 *
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

   const DilutionSchemeInfo *dilPtr;
   const QuarkInfo *qactionPtr;
   const InverterInfo *invertPtr;


       // storage and/or references to internal data

   vector<std::string> fileNames;            // files to handle
   int fileMode;                // 0 = read_only, 1 = must_exist, 
                                // 2 = create_if_not_there, 3 = overwrite
   map<Key,int> fileMap;                    // key -> file index
   map<Key,SerialDBData<LatticeFermion>* > m_storage; // storage of source/sinks


       // Prevent copying ... handler might contain large
       // amounts of data

   QuarkSourceSinkHandler(const QuarkSourceSinkHandler&);
   QuarkSourceSinkHandler& operator=(const QuarkSourceSinkHandler&);



 public:


   QuarkSourceSinkHandler();

   QuarkSourceSinkHandler(XMLReader& xml_in);

   void setInfo(XMLReader& xml_in);

   ~QuarkSourceSinkHandler();

   void clear();

   bool isInfoSet() const;

   string outputInfo() const;

   string getHeader() const;



   const GaugeConfigurationInfo& getGaugeConfigurationInfo() const;

   const FieldSmearingInfo& getFieldSmearingInfo() const;

   const DilutionSchemeInfo& getDilutionSchemeInfo() const; 

   const QuarkInfo& getQuarkInfo() const;

   const InverterInfo& getInverterInfo() const;

   int getNumberOfDilutionProjectors() const;



        // compute for all dilution indices and dump out to file
        // specified by "file_index";  compute sources for all source
        // times, sinks for just one source time but all sink times
        // (file_index == -1 means use the **last** file in the list)

   void computeSource(XMLReader& xml_in, int file_index = -1);

   void computeSource(const LaphNoiseInfo& noise, int file_index = -1);

   void computeSink(XMLReader& xml_in, int file_index = -1);

   void computeSink(const LaphNoiseInfo& noise, int source_time, int file_index = -1);


//   void eraseSink(XMLReader& xml_in);
   
//   void eraseSink(const LaphNoiseInfo& noise, int source_time);
   
//   void eraseSource(XMLReader& xml_in);
   
//   void eraseSource(const LaphNoiseInfo& noise);
   
   

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
   void destroy();
   void set_info(XMLReader& xml_info);
   void set_file_names(XMLReader& record_xml);
   void set_info_helper(XMLReader& xml_info);
   void set_info_helper(const std::string& header);
   void set_info_from_file(XMLReader& xml_in);
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
