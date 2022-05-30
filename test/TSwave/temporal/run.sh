#!/bin/bash
LNS3D_DIR=../../..
FSC_DIR=../../../../fsc
#
# Make the mean profile
#
$FSC_DIR/fsc < fsc.inp
tail -n +2 cprofile.dat > profile.0
tail -n +2 cfirst.dat > first.0
tail -n +2 csecond.dat > second.0
#
# Make the mesh
#
$LNS3D_DIR/util/genmesh < genmesh.inp
#
# Make mean flow from FSC profile
#
$LNS3D_DIR/util/mkmean < mkmean.inp
#
# Make initial condition from time.1 eigenfunction
#
$LNS3D_DIR/util/mkeig3d < mkeig3d.inp
#
# Run the simulation
#
$LNS3D_DIR/src/lns3d < lns3d.inp
#
# Postprocess
#
./mpost.sh output.R.*
#
exit 0
