#!/bin/csh -f
foreach i ($argv[1-])
  echo Processing $i
  ../../util/npost $i 
end
