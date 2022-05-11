#!/bin/bash
\rm -f output.* && \
\rm -f *.plt && \
\rm -f fort.* && \
\rm -f *.dat && \
\rm -f grid.xyz grid3d.xyz out3d.q
exit $?
