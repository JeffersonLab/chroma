#ifndef __inline_psibarpsi_h__
#define __inline_psibarpsi_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma
{
  namespace InlinePsiBarPsiEnv
  {
    extern const std::string name;
	bool registerAll();
  
  
  //! Parameter structure
  /*! \ingroup inlinepbp */
  struct Params
  {
    Params();
	Params(XMLReader& xml_in, const std::string& path);
	void write(XMLWriter& xml_out, const std::string& path);
	
	unsigned long frequency;
	
	struct Param_t
	{
	  GroupXML_t		fermact;
	  GroupXML_t		invParam;
	  
	  multi1d<Real>		mass;
	  int				ichiral;
	} param;
	
	struct NamedObject_t
	{
	  std::string		gauge_id;
	} named_obj;
  };
  
  class InlineMeas : public AbsInlineMeasurement
  {
  public:
    ~InlineMeas() {}
	InlineMeas(const Params& p) : params(p) {}
	InlineMeas(const InlineMeas& p) : params(p.params) {}
	
	unsigned long getFrequency(void) const { return params.frequency; }
	
	void operator()(unsigned long update_no,
			XMLWriter& xml_out);
  protected:
    void func(unsigned long update_no,
			XMLWriter& xml_out);

  private:
    Params params; 
  };

} // end namespace

};

#endif

  

