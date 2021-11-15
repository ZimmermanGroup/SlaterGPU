#!/bin/bash

for d in {1..5}
do
    cd geom_$d
    ../test_geom_timings > ../timings_${d}.log
    cd ..
done