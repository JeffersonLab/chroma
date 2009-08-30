// $Id: inline_stoch_laph_quark_w.cc,v 3.12 2009-08-30 23:56:21 colin Exp $
/*! \file
 * \brief Compute the laph-diluted quark sources and sinks. Write them 
 *  out to db files.  Uses a QuarkSourceSinkHandler.
 *
 * Propagator calculation on laph diluted sources
 */


#include "inline_stoch_laph_quark_w.h"
#include "chroma.h"

namespace Chroma {
  using namespace LaphEnv;
  namespace InlineStochLaphQuarkEnv {

    //  The crucial create measurement routine. Must be in the *.cc
    //  so that it is local to this file.  Dynamically allocates
    //  and instantiates an object of our class "StochLaphQuarkInlineMeas".

AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
                                        const std::string& path) 
{
 return new StochLaphQuarkInlineMeas(xml_in, path);
}

    //  The name of this inline measurement.   This is the name 
    //  with which the createMeasurement function is associated in the 
    //  Object Factory. You must include this name in the XML input 
    //  to Chroma through
    //     <InlineMeasurements>
    //        <elem>
    //            <Name> STOCH_LAPH_QUARK </Name>
    //             ...
    //        </elem>
    //    </InlineMeasurements>

const std::string name = "STOCH_LAPH_QUARK";

    // Registration boolean hidden in anonymous namespace.
namespace {
   bool registered = false;
}

    // Register all the factories.  This function may be called many
    // times by other measurements, so we only want to register this
    // inline measurement once.  Hence, the use of the "registered"
    // boolean above (which must be hidden in an anonymous namespace).

bool registerAll() 
{
 bool success = true; 
 if (!registered){
    success &= TheInlineMeasurementFactory::Instance().registerObject(
                      name, createMeasurement);
    registered = true;}
 return success;
}

// *********************************************************************
	
     // XML input must have form:
     //
     //   <QuarkSourceSinkInfo> ... </QuarkSourceSinkInfo>

     //   <LaphNoiseList> ... </LaphNoiseList>
     //   <SourceTimeList> ... </SourceTimeList> 

     // Inside the <LaphNoiseList> should be one or more <LaphNoise>
     // tags.
     // Inside the <SourceTimeList> should be either
     //     <Values> 2 5 8 </Values>
     // or  <All> </All>





void StochLaphQuarkInlineMeas::clearSinkComputations()
{
 sinkComputations.clear();
}

void StochLaphQuarkInlineMeas::clearSourceComputations()
{
 sourceComputations.clear();
}
      
void StochLaphQuarkInlineMeas::setSinkComputations(int TimeExtent)
{
 if (!sinkComputations.empty()) sinkComputations.clear();

 if (xml_tag_count(xml_rdr,"SinkComputations")==1){
    XMLReader xmlrd(xml_rdr,"./descendant-or-self::SinkComputations");

    if (xml_tag_count(xmlrd,"NoiseList_TimeList_OneFile")==1){
       XMLReader xmlr(xmlrd,"./descendant-or-self::NoiseList_TimeList_OneFile");
       multi1d<int> source_times;
       if (xml_tag_count(xmlr,"SourceTimeList")==1){
          if (xml_tag_count(xmlr,"SourceTimeList/All")==1){
             source_times.resize(TimeExtent);
             for (int t=0;t<TimeExtent;t++) source_times[t]=t;}
          else
             xmlread(xmlr,"SourceTimeList/Values",source_times,
                     "STOCH_LAPH_QUARK");}
       int num_noises=xml_tag_count(xmlr,"LaphNoiseList/LaphNoiseInfo");
       for (int k=1;k<=num_noises;k++){
          ostringstream path;
          path << "./descendant::LaphNoiseList/LaphNoiseInfo["<<k<<"]";
          XMLReader xml_noise(xmlr,path.str());
          LaphNoiseInfo aNoise(xml_noise);
          for (int t=0;t<source_times.size();t++){
             sinkComputations.push_back(
                  SinkComputation(aNoise,source_times[t]));}}}

    if (xml_tag_count(xmlrd,"ComputationList")==1){
       XMLReader xmlr(xmlrd,"./descendant-or-self::ComputationList");
       int ncomputations=xml_tag_count(xmlr,"Computation");
       for (int k=1;k<=ncomputations;k++){
          ostringstream path;
          path << "./descendant::Computation["<<k<<"]";
          XMLReader xml_comp(xmlr,path.str());
          LaphNoiseInfo aNoise(xml_comp);
          int source_time;
          xmlread(xml_comp,"SourceTime",source_time,"STOCH_LAPH_QUARK");
          sinkComputations.push_back(
                SinkComputation(aNoise,source_time));}}

    }

 QDPIO::cout << endl << "STOCH_LAPH_QUARK sink computations:"<<endl;
 QDPIO::cout << " Number of sink computations = "<<sinkComputations.size()<<endl;
 int count=0;
 for (list<SinkComputation>::const_iterator it=sinkComputations.begin();
      it!=sinkComputations.end();count++,it++){
    QDPIO::cout <<endl<< "SinkComputation "<<count<<":"<<endl;
    QDPIO::cout << it->Noise.output();
    QDPIO::cout << "<SourceTime>"<<it->SourceTime<<"</SourceTime>"<<endl;}

}


void StochLaphQuarkInlineMeas::setSourceComputations()
{
 if (!sourceComputations.empty()) sourceComputations.clear();

 if (xml_tag_count(xml_rdr,"SourceComputations")==1){
    XMLReader xmlrd(xml_rdr,"./descendant-or-self::SourceComputations");

    if (xml_tag_count(xmlrd,"ComputationList")==1){
       XMLReader xmlr(xmlrd,"./descendant-or-self::ComputationList");
       int ncomputations=xml_tag_count(xmlr,"Computation");
       for (int k=1;k<=ncomputations;k++){
          ostringstream path;
          path << "./descendant::Computation["<<k<<"]";
          XMLReader xml_comp(xmlr,path.str());
          LaphNoiseInfo aNoise(xml_comp);
          sourceComputations.push_back(
                SourceComputation(aNoise));}}
    }

 QDPIO::cout << endl << "STOCH_LAPH_QUARK source computations:"<<endl;
 QDPIO::cout << " Number of source computations = "<<sourceComputations.size()<<endl;
 int count=0;
 for (list<SourceComputation>::const_iterator it=sourceComputations.begin();
      it!=sourceComputations.end();count++,it++){
    QDPIO::cout <<endl<< "SourceComputation "<<count<<":"<<endl;
    QDPIO::cout << it->Noise.output();}

}


// *********************************************************************
	
     // Subroutine which does all of the work!!  Input parameters
     // must be as shown (specified by Chroma).  Actual input to
     // this routine is through the private data member
     //     XMLReader xlm_rdr


void StochLaphQuarkInlineMeas::operator()(unsigned long update_no,
                                          XMLWriter& xml_out) 
{

    // create the handler and set up the info from the
    // XML <QuarkSourceSinkInfo> tag

 QuarkSourceSinkHandler Q(xml_rdr);

 QDPIO::cout << endl <<"Info initialized in QuarkSourceSinkHandler"<<endl;
 {XMLBufferWriter xmlout;
 Q.getFileMap(xmlout);
 cout << xmlout.str()<<endl;}

    // set or compute the Laph eigenvectors (smears the gauge field as needed)

 Q.setLaphEigenvectors();

    // read the list of computations (noises, time sources, file indices)
    // from xml_rdr and store in the "Computations" data member

 setSourceComputations();
 setSinkComputations(Q.getGaugeConfigurationInfo().getTimeExtent());

 START_CODE();
 StopWatch outer;
 outer.start();

 int count=0;
 for (list<SourceComputation>::const_iterator it=sourceComputations.begin();
      it!=sourceComputations.end();count++,it++){
    QDPIO::cout <<endl<< "Now starting source computation "<<count<<":"<<endl;
    Q.computeSource(it->Noise);}

 count=0;
 for (list<SinkComputation>::const_iterator it=sinkComputations.begin();
      it!=sinkComputations.end();count++,it++){
    QDPIO::cout <<endl<< "Now starting sink computation "<<count<<":"<<endl;
    Q.computeSink(it->Noise,it->SourceTime);}

 outer.stop();
 QDPIO::cout << name << ": total time = " << outer.getTimeInSeconds() 
             << " secs" << endl;
 QDPIO::cout << name << ": ran successfully" << endl;

 END_CODE();
} 

// ******************************************************************
  }
} // namespace Chroma
