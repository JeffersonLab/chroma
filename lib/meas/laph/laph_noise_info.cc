#include "laph_noise_info.h"
#include "xml_help.h"
using namespace std;

namespace Chroma {
  namespace LaphEnv {


// *************************************************************

   // XMLReader constructor

LaphNoiseInfo::LaphNoiseInfo(XMLReader& xml_in)
{
 if (xml_tag_count(xml_in,"LaphNoiseInfo")!=1){
    QDPIO::cerr << "Bad XML input to LaphNoiseInfo"<<endl;
    QDPIO::cerr << "Expected one <LaphNoiseInfo> tag"<<endl;
    QDP_abort(1);}
 XMLReader xmlr(xml_in, "./descendant-or-self::LaphNoiseInfo");
 extract_info_from_reader(xmlr);
}

void LaphNoiseInfo::extract_info_from_reader(XMLReader& xml_in)
{
 xmlread(xml_in,"ZNGroup", znGroup, "LaphNoiseInfo");
 xmlread(xml_in,"seed0", s0, "LaphNoiseInfo");
 xmlread(xml_in,"seed1", s1, "LaphNoiseInfo");
 xmlread(xml_in,"seed2", s2, "LaphNoiseInfo");
 xmlread(xml_in,"seed3", s3, "LaphNoiseInfo");
 check_assignment();
 currTraj=1000;
 calc_seed();
}

void LaphNoiseInfo::check_assignment()
{
 if  ((znGroup<2)||(s0<0)||(s0>4095)||(s1<0)||(s1>4095)
    ||(s2<0)||(s2>4095)||(s3<0)||(s3>2047)){
    QDPIO::cerr << "improper initialization of LaphNoiseInfo"<<endl;
    QDPIO::cerr << "seed0,seed1,seed2 must have value 0 to 4095"<<endl;
    QDPIO::cerr << "seed3 must have value 0 to 2047"<<endl;
    QDPIO::cerr << "ZNGroup must have integer value >= 2"<<endl;
    QDP_abort(1);}
}

 // *************************************************************

   // This version of the constructor assumes that header information
   // from a quark_source_sink file, for example, is passed in.

LaphNoiseInfo::LaphNoiseInfo(const string& header)
{
 string noise_header;
 extract_xml_element(header,"LaphNoiseInfo",noise_header,
                     "LaphNoiseInfo");
 stringstream tmp;
 tmp << noise_header;
 XMLReader xmlr0(tmp);
 XMLReader xmlr(xmlr0,"/LaphNoiseInfo");  
 extract_info_from_reader(xmlr);
}



  // ************************************************************

    // copy constructor

LaphNoiseInfo::LaphNoiseInfo(const LaphNoiseInfo& in)
              : s0(in.s0), s1(in.s1), s2(in.s2), s3(in.s3),
                znGroup(in.znGroup), currTraj(in.currTraj) 
{
 rngseed=in.rngseed;  // no copy constructor for this type!!
}


LaphNoiseInfo& LaphNoiseInfo::operator=(const LaphNoiseInfo& in)
{
 s0=in.s0; s1=in.s1; s2=in.s2; s3=in.s3;
 znGroup=in.znGroup;
 currTraj=in.currTraj;
 rngseed=in.rngseed;
 return *this;
}

void LaphNoiseInfo::checkEqual(const LaphNoiseInfo& in) const
{
 if  ((s0!=in.s0)||(s1!=in.s1)||(s2!=in.s2)||(s3!=in.s3)
    ||(znGroup!=in.znGroup))
    throw string("LaphNoiseInfo does not checkEqual");
}

bool LaphNoiseInfo::operator==(const LaphNoiseInfo& in) const
{
 return ((s0==in.s0)&&(s1==in.s1)&&(s2==in.s2)&&(s3==in.s3)
        &&(znGroup==in.znGroup));
}

bool LaphNoiseInfo::operator<(const LaphNoiseInfo& in) const
{
 return     ((s0<in.s0) || ((s0==in.s0) 
         && ((s1<in.s1) || ((s1==in.s1) 
         && ((s2<in.s2) || ((s2==in.s2) 
         && ((s3<in.s3) || ((s3==in.s3) 
         &&  (znGroup<in.znGroup)))))))));
}


// **********************************************************


Seed LaphNoiseInfo::getSeed(const GaugeConfigurationInfo& G) const
{
 int traj_num=G.getTrajNum();
 if (traj_num==currTraj) return rngseed;
 currTraj=traj_num;
 calc_seed();
 return rngseed;
}

  // perform linear congruential on the four integers
  // i0,i1,i2 are mod 2^12, and i3 is mod 2^11
  // multiplier,shift,mod chosen to satisfy
  //    gcd of shift and mod is 1
  //    multiplier-1 is divisible by all prime factors of mod
  //    multiplier-1 is multiple of 4 when mod is multiple of 4
  //    multiplier*mod+shift must not overflow 2^31

void LaphNoiseInfo::lcongr(int& i0, int& i1, int& i2, int& i3) const
{
 i0=(16385*i0+13) % 4096;   // mod of power of 2 is bit shift...fast
 i1=(16401*i1+27) % 4096;
 i2=(16417*i2+11) % 4096;
 i3=(16449*i3+47) % 2048;
}


void LaphNoiseInfo::calc_seed() const
{
 int nHits=1;
 int trajOffset=1000;
 int i0=s0, i1=s1, i2=s2, i3=s3;

   // run four sequences of linear congruential RNG
   // based on the trajectory number, the offset, and nHits
   //   -> if currTraj>trajOffset, do nHits as traj number increases by 1
   //   -> if currTraj<trajOffset, do nHits as traj number decreases by 1
   //        but jumble the four integers differently
   // RNG used keep i0,i1,i2 twelve bits or below, and i3 eleven bits or below

 for (int k=0;k<nHits*(trajOffset-currTraj);k++)
    lcongr(i1,i2,i0,i3);
 for (int k=0;k<nHits*(currTraj-trajOffset);k++)
    lcongr(i0,i1,i2,i3);


     // put results into the 47 bit seed

 rngseed = i3;
 rngseed = (rngseed<<12) | i2;
 rngseed = (rngseed<<12) | i1;
 rngseed = (rngseed<<12) | i0;

}


string LaphNoiseInfo::output(int indent) const
{
 string pad(3*indent,' ');
 ostringstream oss;
 oss << pad << "<LaphNoiseInfo>"<<endl;
 oss << pad << "  <ZNGroup> " << znGroup << " </ZNGroup>"<<endl;
 oss << pad << "  <seed0> " << s0 << " </seed0>"<<endl;
 oss << pad << "  <seed1> " << s1 << " </seed1>"<<endl;
 oss << pad << "  <seed2> " << s2 << " </seed2>"<<endl;
 oss << pad << "  <seed3> " << s3 << " </seed3>"<<endl;
 oss << pad << "</LaphNoiseInfo>"<<endl;
 return oss.str();
}

void LaphNoiseInfo::output(XMLWriter& xmlout) const
{
 push(xmlout,"LaphNoiseInfo");
 write(xmlout,"ZNGroup",znGroup);
 write(xmlout,"seed0",s0);
 write(xmlout,"seed1",s1);
 write(xmlout,"seed2",s2);
 write(xmlout,"seed3",s3);
 pop(xmlout);
}

void LaphNoiseInfo::binaryWrite(BinaryWriter& out) const
{
 try{
    write(out,znGroup);
    write(out,s0);
    write(out,s1);
    write(out,s2);
    write(out,s3);}
 catch(...){
    QDPIO::cerr << "failed to binary write LaphNoiseInfo"<<endl;
    QDP_abort(1);}
}


void LaphNoiseInfo::binaryRead(BinaryReader& in)
{
 try{
    read(in,znGroup);
    read(in,s0);
    read(in,s1);
    read(in,s2);
    read(in,s3);}
 catch(...){
    QDPIO::cerr << "failed to binary read LaphNoiseInfo"<<endl;
    QDP_abort(1);}
 check_assignment();
}

// *************************************************************
  }
}
