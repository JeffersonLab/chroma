#include "laph_noise_info.h"
#include "xml_help.h"
using namespace std;

namespace Chroma {
  namespace LaphEnv {



   // constructor gets input from XMLReader

LaphNoiseInfo::LaphNoiseInfo(XMLReader& xml_rdr)
{
 try{
    XMLReader xml_in(xml_rdr, "//laph_noise"); 
    read(xml_in,"./seed0", s0 );
    read(xml_in,"./seed1", s1 );
    read(xml_in,"./seed2", s2 );
    read(xml_in,"./seed3", s3 );
    read(xml_in,"./traj_offset", trajOffset );
    read(xml_in,"./nhits", nHits );
    }
 catch(const string& err){
    QDPIO::cerr << "could not initialize LaphNoiseInfo from XML input"<<endl;
    QDP_abort(1);}
 if (nHits<1){
    QDPIO::cerr << "nHits > 0 required in LaphNoiseInfo"<<endl;
    QDP_abort(1);}
 
 currTraj=0;
 rngseed = s0;
 /*
		rngseed <<= 12;
 rngseed |= s1;
 rngseed <<= 12;
 rngseed |= s2;
 rngseed <<= 12;
 rngseed |= s3; 
*/
 }

Seed LaphNoiseInfo::getSeed(const GaugeConfigurationInfo& G) const
{
 int traj_num=G.getTrajNum();
 if (traj_num==currTraj) return rngseed;
 unsigned int ss0=s0, ss1=s1, ss2=s2, ss3=s3;
 for (int i=0;i<(nHits+1)*(trajOffset-traj_num);i++)
    lcongr(ss0,ss1,ss2,ss3);
 for (int i=0;i<nHits*(traj_num-trajOffset);i++)
    lcongr(ss0,ss1,ss2,ss3);
 currTraj=traj_num;
 rngseed = 0;
 rngseed = ss0;
/* 
 rngseed <<= 12;
 rngseed |= ss1;
 rngseed <<= 12;
 rngseed |= ss2;
 rngseed <<= 12;
 rngseed |= ss3; 
 */
 return rngseed;
}


void LaphNoiseInfo::lcongr(unsigned int& ss0, unsigned int& ss1,
                           unsigned int& ss2, unsigned int& ss3) const
{
 ss0=(ran_mult*ss0+ran_shift) % ran_mod;
 ss1=(ran_mult*ss1+ran_shift) % ran_mod;
 ss2=(ran_mult*ss2+ran_shift) % ran_mod;
 ss3=(ran_mult*ss3+ran_shift) % ran_mod;
}


string LaphNoiseInfo::output() const
{
 ostringstream oss;
 oss << "<laph_noise>"<<endl;
 oss << "  <seed0> " << s0 << " </seed0>"<<endl;
 oss << "  <seed1> " << s1 << " </seed1>"<<endl;
 oss << "  <seed2> " << s2 << " </seed2>"<<endl;
 oss << "  <seed3> " << s3 << " </seed3>"<<endl;
 oss << "  <traj_offset> " << trajOffset << " </traj_offset>"<<endl;
 oss << "  <nhits> " << nHits << " </nhits>"<<endl;
 oss << "</laph_noise>"<<endl;
 return oss.str();
}


// *************************************************************
  }
}
