#!/bin/csh -f
set echo on
conv -g grid.dat
mv grid.bin tmp.g
foreach i ($argv[1-])
  ~/codes/util/lpost3d -ij $i > /dev/null
  set j = $i:r
  conv $j:r.q.$i:e
  cp $j:r.b.$i:e tmp.q
  preplot tmp prim.$i:e.plt -plot3d -b -3dw > /dev/null
end
rm tmp.g tmp.q
