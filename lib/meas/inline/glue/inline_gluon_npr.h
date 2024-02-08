// -*- C++ -*-
/*! \file
 *  \brief Inline naive topological charge
 *
 * Author: Joe Karpie
 */

#ifndef INLINE_GLUON_NPR_H
#define INLINE_GLUON_NPR_H

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

#include <cmath>

namespace Chroma 
{ 
  /*! \ingroup inlineglue */
  namespace InlineGluonNPREnv 
  {
    bool registerAll();

    /*! \ingroup inlineglue */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);

      unsigned long frequency;

      struct Param_t
      {
	Real          r1;
	Real          r2;
        int         rho;
	multi1d<int>  p;
      } param;

      struct NamedObject_t
      {
	std::string   gauge_id;
	std::string   obs_id;
      } named_obj;
    };
    

    /*! \ingroup inlineglue */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      void operator()(unsigned long update_no,
		      XMLWriter& xml_out); 

    protected:
      void func(const unsigned long update_no, 
		XMLWriter& xml_out);

    private:
      Params params;
      Set set_r1;
      Set set_r2;
    };

    class CloseRangeFunc  : public SetFunc
    {
    public:
      CloseRangeFunc(Real rr, multi1d<int> xx) { r2 = rr*rr; x0 = xx; }
      int operator() (const multi1d<int>& coordinate) const 
      {
        multi1d<int> ndim=Layout::lattSize();
        Real x2 = 0;
        for(int i = 0; i < Nd ; i++)
        {
          int x = abs(coordinate[i] - x0[i] );
          x = fmin(x,(ndim[i] - x));
          x2 += x*x;
        }
        if(toBool(r2 >= x2))
          return 0;
        else
          return 1;
      };
      int numSubsets() const {return 2;} ;
      
    private:
      CloseRangeFunc() {}  // hide default constructor

      Real r2;
      multi1d<int> x0;
    };


  }

}

#endif
