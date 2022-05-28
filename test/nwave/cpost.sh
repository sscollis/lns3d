#!/bin/csh -f
set echo on
conv -g grid.dat
cp grid.bin tmp.g
foreach i ($argv[1-])
  npost -rc -D5 $i > /dev/null
  set j = $i:r
  conv $j:r.q.$i:e
  conv $j:r.d.$i:e
  cp $j:r.b.$i:e tmp.q
  preplot tmp prim.$i:e.plt -plot3d -b -3dw > /dev/null
  cp $j:r.d.$i:e.bin tmp.q
  preplot tmp deriv.$i:e.plt -plot3d -b -3dw > /dev/null
  rm $j:r.q.$i:e $j:r.d.$i:e $j:r.b.$i:e $j:r.d.$i:e.bin
end
rm tmp.g tmp.q
