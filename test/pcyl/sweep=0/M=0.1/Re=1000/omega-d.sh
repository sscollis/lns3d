#!/bin/bash
awk 'NR>2 {print sqrt(2*$1*1000), -$3*(2*$1+1)/sqrt(2*$1*1000)}' wall.dat > omega-d.dat
exit 0
