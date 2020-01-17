#!/bin/bash

pno=$1

(cd mercury_$pno/; ./mercury6_2_$pno) &
while pgrep -x mercury6_2_$pno > /dev/null
do
    if grep -e was\ hit\ by -e ejected\ at -e collided\ with\ the\ central\ body mercury_$pno/info.out
    then
        kill `pidof mercury6_2_$pno`
    fi
done

(cd mercury_$pno/; ./element6_2)
