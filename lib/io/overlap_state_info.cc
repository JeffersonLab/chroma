
#include "chromabase.h"
#include "overlap_state_info.h"

namespace Chroma { 

OverlapStateInfo::OverlapStateInfo(void) 
{
  initedP = false;
  ApproxMin = 0;
  ApproxMax = 0; 
  NWilsVec = 0;
  load_eigenP = false;
 
  
  eigen_io.eigen_file="";
  eigen_io.eigen_volfmt = QDPIO_SINGLEFILE;

  ritzery.Neig = 0;
  ritzery.RsdR = 0;
  ritzery.RsdA = 0;
  ritzery.RsdZero = 0;
  ritzery.ProjApsiP = false;
  ritzery.Ndummy = 0;
  ritzery.GammaFactor = 0;
  ritzery.MaxKS =0;
  ritzery.MinKSIter =0;
  ritzery.MaxKSIter =0;
  ritzery.MaxCG = 0;
  ritzery.Nrenorm = 0;
  
}

void read(XMLReader& xml_in, const string& path, OverlapStateInfo& info)
{
  bool initedP;
  Real ApproxMin;
  Real ApproxMax;
  int  NWilsVec;
  bool load_eigenP;
  EigenIO_t eigen_io;
  RitzParams_t ritzery;

  initedP = false;
  ApproxMin = 0;
  ApproxMax = 0; 
  NWilsVec = 0;
  load_eigenP = false;
  
  eigen_io.eigen_file="";
  eigen_io.eigen_volfmt = QDPIO_SINGLEFILE;

  ritzery.Neig = 0;
  ritzery.RsdR = 0;
  ritzery.RsdA = 0;
  ritzery.RsdZero = 0;
  ritzery.ProjApsiP = false;
  ritzery.Ndummy = 0;
  ritzery.GammaFactor = 0;
  ritzery.MaxKS =0;
  ritzery.MinKSIter =0;
  ritzery.MaxKSIter =0;
  ritzery.MaxCG = 0;
  ritzery.Nrenorm = 0;
  
  XMLReader reader(xml_in, path);
  try 
  {
    if( reader.count("NWilsVec") != 0 )
      read(reader, "NWilsVec", NWilsVec);
  }
  catch(const string& e) { 
    QDPIO::cerr << "Caught exception : " << e << endl;
    QDP_abort(1);
  }


  if( NWilsVec == 0 ) 
  {
    // No eigenbeasties
    eigen_io.eigen_file="";
    eigen_io.eigen_volfmt = QDPIO_SINGLEFILE;
    
    // Try reading approx min and optional approx max
    // if approx max is not specified set it to 2*Nd
    try { 
      if( reader.count("ApproxMin") == 0 )
	ApproxMin = 0.0;
      else
	read(reader, "ApproxMin", ApproxMin);

      if( reader.count("ApproxMax") == 0 )
	ApproxMax = 2*Nd;
      else
	read(reader, "ApproxMax", ApproxMax);
    }
    catch( const string& e) {
      QDPIO::cerr << "Caught exception : " << e << endl;
    }
  }
  else {
    // We have eigenbeasties
    ApproxMin =0;
    ApproxMax =0;

    if( reader.count("Eig") == 1 ) {
      if( reader.count("Ritz") == 1 ) {
	QDPIO::cerr << "Cannot specify both Eig and Ritz " << endl;
	QDP_abort(1);
      }


      // Read in the eigenvector IO params
      try { 
	read(reader, "Eig", eigen_io);
      }
      catch( const string& e ) { 
	QDPIO::cerr << "Caught exception: " << e << endl;
      }
      load_eigenP = true;

    }
    else if ( reader.count("Ritz") == 1 ) {
      if( reader.count("Eig") == 1 ) { 
	QDPIO::cerr << "Cannot specify both Eig and Ritz " << endl;
	QDP_abort(1);
      }
      
      try {
	XMLReader ritzreader(reader, "Ritz");

	read(ritzreader, "RsdR", ritzery.RsdR);
	read(ritzreader, "RsdA", ritzery.RsdA);
	read(ritzreader, "RsdZero", ritzery.RsdZero);
	read(ritzreader, "ProjApsiP", ritzery.ProjApsiP);
	read(ritzreader, "Ndummy",    ritzery.Ndummy);
	read(ritzreader, "GammaFactor", ritzery.GammaFactor);
	read(ritzreader, "MaxKS", ritzery.MaxKS);
        read(ritzreader, "MaxCG", ritzery.MaxCG);
	read(ritzreader, "MinKSIter", ritzery.MinKSIter);
	read(ritzreader, "MaxKSIter", ritzery.MaxKSIter);
	read(ritzreader, "Nrenorm", ritzery.Nrenorm);
      }    
      catch( const string& e ) { 
	QDPIO::cerr << "Caught exception: " << e << endl;
	QDP_abort(1);
      }

      ritzery.Neig = NWilsVec;
      load_eigenP = false;
    }
    else { 
      QDPIO::cerr << "Must specify either Eig for loadable eigenvalues or "
		  << "Ritz Parameters for compuing eigenvalues" << endl;
      QDP_abort(1);
    }

  }

  info.init(ApproxMin,
	    ApproxMax,
	    NWilsVec,
	    load_eigenP,
	    eigen_io,
	    ritzery);

}

void write(XMLWriter& xml_out, const string& path, const OverlapStateInfo& info)
{

  if( path != "." ) { 
    push(xml_out, path);
  }

  write( xml_out, "NWilsVec", info.getNWilsVec());

  if( info.getNWilsVec() == 0 ) {    
    write(xml_out, "ApproxMin", info.getApproxMin());
    write(xml_out, "ApproxMax", info.getApproxMax());
  }
  else {
    if( info.loadEigVec() ) { 
      write(xml_out, "Eig", info.getEigenIO());
    }
    else if ( info.computeEigVec() ) { 
      
      const RitzParams_t& ritzery = info.getRitzParams();

      push(xml_out, "Ritz");
      write(xml_out, "RsdR", ritzery.RsdR);
      write(xml_out, "RsdA", ritzery.RsdA);
      write(xml_out, "RsdZero", ritzery.RsdZero);
      write(xml_out, "ProjApsiP", ritzery.ProjApsiP);
      write(xml_out, "Ndummy",    ritzery.Ndummy);
      write(xml_out, "GammaFactor", ritzery.GammaFactor);
      write(xml_out, "MaxKS", ritzery.MaxKS);
      write(xml_out, "MaxCG", ritzery.MaxCG);
      write(xml_out, "MinKSIter", ritzery.MinKSIter);
      write(xml_out, "MaxKSIter", ritzery.MaxKSIter);
      write(xml_out, "Nrenorm", ritzery.Nrenorm);
      pop(xml_out);

    }
    else { 
      QDPIO::cerr << "Must specify either Eig for loadable eigenvalues or "
		  << "Ritz Parameters for compuing eigenvalues" << endl;
    }

  }

  if(path != "." ) { 
    pop(xml_out);
  }

}

};
