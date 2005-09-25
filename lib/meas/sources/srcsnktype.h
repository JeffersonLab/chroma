// -*- C++ -*-
// $Id: srcsnktype.h,v 2.0 2005-09-25 21:04:40 edwards Exp $
/*! \file
 *  \brief Source and sink types
 */

#ifndef __srcsnktype_h__
#define __srcsnktype_h__

namespace Chroma {
  enum SourceType
    {
      SRC_TYPE_POINT_SOURCE,
      SRC_TYPE_WALL_SOURCE,
      SRC_TYPE_SHELL_SOURCE,
      SRC_TYPE_BNDST_SOURCE,
      SRC_TYPE_POINT_AND_BNDST_SOURCE,
      SRC_TYPE_SHELL_AND_BNDST_SOURCE,
      SRC_TYPE_POINT_AND_SHELL_AND_BNDST_SOURCE,
      SRC_TYPE_RAND_Z2_WALL_SOURCE,
      SRC_TYPE_RAND_U1_WALL_SOURCE
    };

  enum SinkType
    {
      SNK_TYPE_POINT_SINK,
      SNK_TYPE_WALL_SINK,
      SNK_TYPE_POINT_AND_WALL_SINK,
      SNK_TYPE_SHELL_SINK,
      SNK_TYPE_POINT_AND_SHELL_SINK,
      SNK_TYPE_BNDST_SINK,
      SNK_TYPE_POINT_AND_BNDST_SINK,
      SNK_TYPE_SHELL_AND_BNDST_SINK,
      SNK_TYPE_POINT_AND_SHELL_AND_BNDST_SINK
    };
}

#endif
