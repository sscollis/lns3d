#!/bin/bash
#
# exit if any command fails
#
set -e
#
#=============================================================================
#     L i n e a r i z e d   N a v i e r - S t o k e s   D r i v e r
#=============================================================================
# Author:     S. Scott Collis
# Copyright:  S. Scott Collis(c)
# Date:       08/09/2022
#
# Purpose:    Run a bump receptivity calculation similar to case in Collis'
#             PhD thesis, case 3 of Table 5.2
#=============================================================================
#
# setup paths (your might need to supply alternatives)
#
LNS3D_DIR=${LNS3D_DIR:=$HOME/git/lns3d}
#
# run the case
#
$LNS3D_DIR/src/lns3d < lns3d.inp
#
# visualize
#
./lpost3d.sh output.R.*
#
exit $?
