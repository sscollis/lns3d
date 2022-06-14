#!/bin/bash
\rm -f output.[Rhq].* *.res && \
\rm -f state.* *.ij *.log && \
\rm -f *.jpeg *.gif && \
\rm -f npot.f.0 npot.q.0 *.xyz && \
\rm -f fort.* \
\rm \
amp.dat \
bl3d.dat \
complex.q.0 \
coord.dat \
cp.dat \
echo.dat \
grid.dat \
ic.dat \
lns.dat \
m1.dat \
m11.dat \
mag.dat \
mean.dat \
mean.q.0 \
n11.dat \
real.q.0 \
restart.dat \
sponge.dat \
sponge2.dat \
stat.dat \
wall.dat \
top.dat
exit $?
