#include "field_smearing_handler.h"
#include "xml_help.h"
#include <typeinfo>

namespace Chroma {
  namespace LaphEnv {


   // constructors

FieldSmearingHandler::FieldSmearingHandler() 
     : uPtr(0), smearPtr(0), laph_eigvecs(0) {}



FieldSmearingHandler::FieldSmearingHandler(XMLReader& smear_xml)
     : uPtr(0), laph_eigvecs(0),smearPtr(0)
{
 setInfo(smear_xml);
}


void FieldSmearingHandler::setInfo(XMLReader& xml_info)
{
 clear();
 create_handlers();
 try{
    uPtr->setInfo(xml_info);
    if (xml_info.count("//laph_eigvecs_namedobj_id")==1)
       read(xml_info,"//laph_eigvecs_namedobj_id",laph_eigvecs_id);
   else
       laph_eigvecs_id="laph_eigvecs_id";
    smearPtr=new FieldSmearingInfo(xml_info);}
 catch(...){
    QDPIO::cerr << "problem in setInfo in FieldSmearingHandler"<<endl;
    QDP_abort(1);}
 laph_eigvecs_id=tidyString(laph_eigvecs_id);
 QDPIO::cout << uPtr->outputInfo() <<endl;
 QDPIO::cout << "FieldSmearingInfo set in FieldSmearingHandler"<<endl;
}

    // destructor

FieldSmearingHandler::~FieldSmearingHandler()
{
 clear();
}


void FieldSmearingHandler::clear()
{
 try {delete smearPtr;} catch(...) {QDP_abort(1);}
 smearPtr=0;
 laph_eigvecs_id.clear();
 destroy_handlers(); 
 clear_data();  
 QDPIO::cout << "FieldSmearingHandler cleared"<<endl;
}


bool FieldSmearingHandler::isInfoSet() const
{
 return ((smearPtr!=0)&&(uPtr->isInfoSet())&&(!laph_eigvecs_id.empty()));
}

string FieldSmearingHandler::outputInfo() const
{
 ostringstream oss;
 oss << smearPtr->output() << endl;
 oss << uPtr->outputInfo() << endl;
 return oss.str();
}
 

void FieldSmearingHandler::computeLaphEigenvectors()
{
 if (!isInfoSet()){
    QDPIO::cerr << " must setInfo in FieldSmearingHandler before"<<endl
                << "   computing Laph Eigenvectors" << endl;
    QDP_abort(1);}
 if (laph_eigvecs!=0){
    QDPIO::cerr << "clear current Laph Eigenvectors before computing"<<endl;
    QDP_abort(1);}
    
 try{
    TheNamedObjMap::Instance().create< SubsetVectors<LatticeColorVector> >(
                 laph_eigvecs_id);
    }
 catch (std::bad_cast){
    QDPIO::cerr << "caught dynamic cast error in computeLaphEigenvectors" << endl;
    QDP_abort(1);
    }
 catch (const string& err){
    QDPIO::cerr << "error in computeLaphEigenvectors: " << err << endl;
    QDP_abort(1);
    }
 laph_eigvecs = &(TheNamedObjMap::Instance().getData<SubsetVectors<LatticeColorVector> >(
                   laph_eigvecs_id));

 int nEigvecs = smearPtr->getNumberOfLaplacianEigenvectors();
 int nTime = uPtr->getInfo().getTimeExtent();
 int tdir = uPtr->getInfo().getTimeDir();

 laph_eigvecs->getEvectors().resize(nEigvecs);
 laph_eigvecs->getEvalues().resize(nEigvecs);
 for (int k=0;k<nEigvecs;k++)
    laph_eigvecs->getEvalues()[k].weights.resize(nTime);

 QDPIO::cout << "Computation of Laplacian eigenvectors commencing"
             << " in FieldSmearingHandler"<<endl;
 QDPIO::cout << "  number of requested eigenvectors = "<<nEigvecs<<endl;
 QDPIO::cout << "            time extent of lattice = "<<nTime<<endl;

 START_CODE();
 StopWatch rolex;
 rolex.start();


// do the computation, put into laph_eigvecs  

  // ***********************DUMMY*****************************
  
  for (int k =0;k<nEigvecs;k++)
     gaussian(laph_eigvecs->getEvectors()[k]);
 
  
  
  
  // ***********************DUMMY*****************************

 rolex.stop();
 QDPIO::cout << "computeLaphEigenvectors: total time = "
             << rolex.getTimeInSeconds() 
             << " secs" << endl;
 QDPIO::cout << "ran successfully" << endl;

   // put metadata into NamedObjMap

 XMLBufferWriter file_xml;
 push(file_xml, "LaplaceEigvectorInfo");
 pop(file_xml);

 XMLBufferWriter record_xml;
 push(record_xml, "LaplaceEigvectorInfo");
 record_xml << uPtr->getInfo().output();
 record_xml << smearPtr->output();
 push(record_xml, "SubsetVectors");
 for(int i=0;i<nEigvecs;i++){
   push(record_xml, "EigenPair");
   write(record_xml, "EigenPairNumber", i); 
   write(record_xml, "EigenValues", laph_eigvecs->getEvalues()[i].weights); 
   pop(record_xml);
 }
 pop(record_xml);
 pop(record_xml);

 // Write metadata to NamedObjMap
 TheNamedObjMap::Instance().get(laph_eigvecs_id).setFileXML(file_xml);
 TheNamedObjMap::Instance().get(laph_eigvecs_id).setRecordXML(record_xml);

 // dump out some xml log information
 QDPIO::cout << "Laplacian eigenvectors computed in FieldSmearingHandler"<<endl;
 QDPIO::cout << record_xml.str();
}



  // read eigenvectors that we previously calculated and stored
  // in the named object map

void FieldSmearingHandler::setLaphEigenvectors()
{
 if (!isInfoSet()){
    QDPIO::cerr << " must setInfo in FieldSmearingHandler before"<<endl
                << "   setting Laph Eigenvectors" << endl;
    QDP_abort(1);}
 if (laph_eigvecs!=0){
    QDPIO::cerr << "warning: setLaphEigenvalues called when already set"<<endl;
    return;}

 XMLReader source_record_xml;
 QDPIO::cout << "Getting Laph eigenvectors from NamedObjMap" << endl;
 if (!TheNamedObjMap::Instance().check(laph_eigvecs_id)){
    QDPIO::cout << "Laph eigenvectors are not currently in NamedObjMap";
    laph_eigvecs=0;
    return;}
 try{
    TheNamedObjMap::Instance().get(laph_eigvecs_id).getRecordXML(source_record_xml);
    check_set(source_record_xml);
    laph_eigvecs = &(TheNamedObjMap::Instance().getData<SubsetVectors<LatticeColorVector> >(laph_eigvecs_id));
    }
 catch (std::bad_cast){
    QDPIO::cerr << ": caught dynamic cast error" << endl;
    QDP_abort(1);
    }
 catch (const string& err){
    QDPIO::cerr << "error in setLaphEigenvectors: " << err << endl;
    QDP_abort(1);
    }
 QDPIO::cout << "LaphEigenvectors successfully set in FieldSmearingHandler" << endl;
}


const multi1d<Real>& FieldSmearingHandler::getLaphEigenvalues(int t) const
{ 
 if (laph_eigvecs==0){
    QDPIO::cerr << "attempt to getLaphEigenvalues before computed/set"<<endl;
    QDP_abort(1);}
 else if ((t<0)||(t>=uPtr->getInfo().getTimeExtent())){
    QDPIO::cerr << "getLaphEigenvalues failed due to invalid time"<<endl;
    QDP_abort(1);}
 return laph_eigvecs->getEvalues()[t].weights;
}

const multi1d<LatticeColorVector>& FieldSmearingHandler::getLaphEigenvectors() const
{
 if (laph_eigvecs==0){
    QDPIO::cerr << "attempt to getLaphEigenvectors before computed/set"<<endl;
    QDP_abort(1);}
 return laph_eigvecs->getEvectors();
}

           // clear pointer to Laph eigenvectors, and if "namedobj_erase"
           // is true, also remove from the NamedObjMap

void FieldSmearingHandler::clearLaphEigenvectors(bool namedobj_erase)
{
 if (namedobj_erase){
    TheNamedObjMap::Instance().erase(laph_eigvecs_id);
    }
 laph_eigvecs = 0;
 QDPIO::cout << "Laph Eigenvectors cleared";
 if (namedobj_erase) QDPIO::cout << " and erased from NamedObjMap";
 QDPIO::cout << "in FieldSmearingHandler"<<endl;
}









void FieldSmearingHandler::computeSmearedGaugeField()
{
 if (!isInfoSet()){
    QDPIO::cerr << " must setInfo in FieldSmearingHandler before"<<endl
                << "  computing smeared gauge field" << endl;
    QDP_abort(1);}

 double rho = smearPtr->getLinkStapleWeight();  cout << "rho = "<<rho<<endl;
 int niters = smearPtr->getNumberOfLinkIterations();
 int tdir = uPtr->getInfo().getTimeDir();
 int ndir = uPtr->getInfo().getNumberOfDirections();

 QDPIO::cout << "Computation of smeared gauge field commencing in FieldSmearingHandler"<<endl;
 QDPIO::cout << smearPtr->output()<<endl;
 
 cout << "tdir = "<<tdir<<endl;
 cout << "ndir = "<<ndir<<endl;
 multi1d<bool> smear_dirs(ndir);
 for (int i=0;i<ndir;i++) smear_dirs[i] = true;
 smear_dirs[tdir]=false;

 multi2d<Real> rhomat(ndir,ndir);
 for(int mu=0; mu < ndir; mu++) 
 for(int nu=0; nu < ndir; nu++){
    if( mu != nu ) rhomat[mu][nu] = Real(rho);
    else rhomat[mu][nu] = 0;}
 for(int mu=0; mu < ndir; mu++)
    if ( ! smear_dirs[mu] ){
       for(int nu=0; nu < ndir; nu++){
          rhomat[mu][nu] = 0;
          rhomat[nu][mu] = 0;}}

 multi1d<LatticeColorMatrix> utemp(ndir);
 usmear.resize(ndir);

 START_CODE();
 StopWatch rolex;
 rolex.start();  cout << "start fuzzing"<<endl;

 if ((niters%2)==0) usmear=uPtr->getData();
 else{
    utemp=uPtr->getData(); cout << "utemp assigned"<<endl;
    Chroma::Stouting::smear_links(utemp,usmear,smear_dirs,rhomat);}

 for (int k=0;k<(niters/2);k++){ cout << "k = "<<k<<endl;
    Chroma::Stouting::smear_links(usmear,utemp,smear_dirs,rhomat);
    Chroma::Stouting::smear_links(utemp,usmear,smear_dirs,rhomat);}

 rolex.stop();
 END_CODE();
 QDPIO::cout << "computeSmearedGaugeField: total time = "
             << rolex.getTimeInSeconds() 
             << " secs" << endl;
 QDPIO::cout << "ran successfully in FieldSmearingHandler" << endl;
}



const multi1d<LatticeColorMatrix>& FieldSmearingHandler::getSmearedGaugeField() const
{
 if (usmear.size()==0){
    QDPIO::cerr << "error in getSmearedGaugeField: unassigned"<<endl;
    QDP_abort(1);}
 QDPIO::cout << "smeared gauge field set in FieldSmearingHandler"<<endl;
 return usmear;
}

        // clear the smeared gauge field from "usmear"

void FieldSmearingHandler::clearGaugeField()
{ 
 usmear.resize(0);
 QDPIO::cout << "smeared gauge field clear in FieldSmearingHandler"<<endl;
}



void FieldSmearingHandler::clear_data()
{
 clearLaphEigenvectors();
 clearGaugeField();
}


void FieldSmearingHandler::create_handlers()
{
 try{
    uPtr=new GaugeConfigurationHandler;}
 catch(...){
    QDPIO::cerr << "allocation problem in FieldSmearingHandler"<<endl;
    QDP_abort(1);}
}

void FieldSmearingHandler::destroy_handlers()
{
 try {delete uPtr;} catch(...) {QDP_abort(1);}
 uPtr=0;
}



void FieldSmearingHandler::check_set(XMLReader& record_xml)
{
 GaugeConfigurationInfo ucheck(record_xml);
 FieldSmearingInfo scheck(record_xml);
 uPtr->getInfo().matchXMLverbatim(ucheck);
 smearPtr->checkEqual(scheck);
}

// *****************************************************************
  }
}
