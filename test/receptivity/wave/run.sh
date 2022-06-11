#!/bin/bash
LNS3D=$HOME/git/lns3d
$LNS3D/util/genmesh < genmesh.inp && \
$LNS3D/util/mkmean < mkmean.inp   && \
$LNS3D/util/mkdist3d < mkdist3d.inp && \
$LNS3D/src/lns3d < lns3d.inp && \
./lpost3d.sh output.R.*
exit $?
