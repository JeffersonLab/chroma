#!/bin/sh

../../scalar/mainprogs/main/ExampleBuildingBlocks.exe 4 4 4 8 \
  ../test_purgaug.cfg1  szin \
  ../propagator/propagator_0 chroma \
  ../seqsource/seqprop_0_PION chroma G5_B_G5 \
  ../seqsource/seqprop_0_PION chroma G5_B_G5 \
  1 3 \
  "pion_u_%c%1d_%c%1d_%c%1d.bb" \
  "pion_d_%c%1d_%c%1d_%c%1d.bb" \
  output_v1.txt output_v1.xml

#  NULL chroma G5_B_G5 \
