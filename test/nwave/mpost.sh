#!/bin/csh -f
foreach i ($argv[1-])
  npost -d $i 
  mv fort.17 f.$i:e
end
