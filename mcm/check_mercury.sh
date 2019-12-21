#!/bin/bash

(cd mercury/; ./mercury6_2) &
while pgrep -x mercury6_2 > /dev/null
do
    if grep -e was\ hit\ by -e ejected\ at -e collided\ with\ the\ central\ body mercury/info.out
    then
        kill `pidof mercury6_2`
    fi
done

(cd mercury/; ./element6_2)
