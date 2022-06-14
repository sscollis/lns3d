#!/bin/bash
\rm -f output.* *.pot *.dat *.res && \
\rm -f state.* *.ij *.log && \
\rm -f *.jpeg *.gif && \
\rm -f npot.f.0 npot.q.0 *.xyz && \
\rm -f fort.* profile.* first.* second.*
exit $?
