#!/bin/bash

pno=$1

(cd mercury_$pno/; ./mercury_$pno.exe) &
while pgrep -x mercury_$pno.exe > /dev/null
do
    if grep -e was\ hit\ by -e ejected\ at -e collided\ with\ the\ central\ body mercury_$pno/info.out
    then
        kill `pgrep mercury_$pno.exe`
    fi
done

(cd mercury_$pno/; ./element_$pno.exe)
