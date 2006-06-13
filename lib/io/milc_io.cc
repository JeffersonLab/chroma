// $Id: milc_io.cc,v 3.1 2006-06-13 18:16:57 bjoo Exp $

/*! \file
 *  \brief MILC gauge format routines
 */

#include "chromabase.h"
#include "io/milc_io.h"
#include <time.h>

namespace Chroma 
{

  //! Initialize header with default values
  MILCGauge_t::MILCGauge_t()
  {
    nrow = Layout::lattSize();

    time_t now = time(NULL);
    {
      char *tmp = ctime(&now);
      int date_size = strlen(tmp);
      char *datetime = new(nothrow) char[date_size+1];
      if( datetime == 0x0 ) { 
	QDP_error_exit("Unable to allocate datetime in qdp_iogauge.cc\n");
      }

      strcpy(datetime,ctime(&now));

      for(int i=0; i < date_size; ++i)
	if ( datetime[i] == '\n' )
	{
	  datetime[i] = '\0';
	  date_size = i;
	  break;
	}   

      date = datetime;
      delete[] datetime;
    }
  }



  //! Source header read
  void read(XMLReader& xml, const string& path, MILCGauge_t& header)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "date", header.date);
    read(paramtop, "nrow", header.nrow);
  }


  //! Source header writer
  void write(XMLWriter& xml, const string& path, const MILCGauge_t& header)
  {
    push(xml, path);

    write(xml, "date", header.date);
    write(xml, "nrow", header.nrow);

    pop(xml);
  }

}  // end namespace Chroma
