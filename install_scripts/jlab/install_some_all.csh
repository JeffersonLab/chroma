#!/bin/tcsh

ssh -n qcdi01   "cd qcd/chroma; ./install_some.csh"
ssh -n qcdi02   "cd qcd/chroma; ./install_some.csh"
