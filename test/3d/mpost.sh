#!/bin/csh -f
foreach file ($argv)
  echo Processing $file
  ../../util/lpost3d -ij $file > /dev/null
end
