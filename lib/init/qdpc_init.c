#include "init/qdpc_init.h"


#include <qop_qdp.h>
#include <qdp_layout.h>

void initializeQDPC(int *argc, char **argv[])
{
  QDP_initialize(argc, argv);
  QDP_profcontrol(0);
  QDP_set_default_layout(QDP_layout_hyper_eo);
  QDP_set_read_group_size(8);
  QDP_set_write_group_size(8);
}

void finalizeQDPC(void)
{
  QDP_finalize();
}
