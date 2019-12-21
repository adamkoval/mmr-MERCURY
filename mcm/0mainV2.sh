#!/bin/bash

res="$1-$2"  	# 1st AND 2nd ARGS USED TO NAME DIRECTORY FOR NOMINAL RESONANCE
K="$3"  		# 3rd ARG IS THE NUMBER OF RUNS TO PERFORM

# TEST IF NOMINAL RESONANCE i
# DIRECTORY EXISTS, IF NOT, CREATE
if ! test -d ../completed/$res/;then
	mkdir ../completed/$res/
fi

if ! test -d 

# COUNT HOW MANY FILES ALREADY EXIST 
# IN e.g. /4-3/ DIR, e.g. 1,2,3 ALREADY
# DONE, SO /4/ IS CREATED
folders="info input planets"
for folder in $folders; do
	COUNT=$(find ../completed/$res/$folder/* -maxdepth 0 | wc -l)
DIR=$((COUNT+1))
mkdir ../completed/$res/$DIR ../completed/$res/$DIR/info/ ../completed/$res/$DIR/input/ ../completed/$res/$DIR/planets/

# MAIN LOOP BEGINS HERE
for k in $( seq 1 $K )
do
	# RANDOMIZE VARIABLES, PASSING 1st
	# & 2nd ARGS TO FUNCTION
	./randomize.sh "$1" "$2"
	wait
	# CLEAN UP OLD FILES 
	# (REDUNDANT NOW MAYBE)
	./cleanup.sh
	wait

	# LAUNCH MERCURY, CHECK FOR UNSTABLE 
	# OUTCOME WHILE IT IS RUNNING AND KILL
	# IF THIS OCCURS
	(cd mercury/; ./mercury6_2) &
	while pgrep -x mercury6_2 > /dev/null
	do
		if grep -e was\ hit\ by -e ejected\ at -e collided\ with\ the\ central\ body mercury/info.out
		then
			kill `pidof mercury6_2`
		fi
	done

	# LAUNCH ELEMENT6 TO CONVERT XV.OUT
	# INTO READABLE FILES
	(cd mercury/; ./element6_2)

	# MAKE NEW FOLDER FOR THE CURRENT RUN
	# TO STORE PLANET FILES AND COPY FILES INTO IT
	if ! test -d ../completed/$res/$DIR/planets/$k/
	theN
		mkdir -p ../completed/$res/$DIR/planets/$k/ && cp mercury/*.aei ../completed/$res/$DIR/planets/$k/

		# REMOVE HEADERS FROM FILES TO MAKE THEM PYTHON-READABLE
		for j in 1 2
		do
			tail -n +5 ../completed/$res/$DIR/planets/$k/planet$j.aei > ../completed/$res/$DIR/planets/$k/planet$j\b.aei
		done
		
		# COPY INFO.OUT AND BIG.IN INTO NEW DIRECTORIES
		cp mercury/info.out ../completed/$res/$DIR/info/$k\-info.out
		cp mercury/big.in ../completed/$res/$DIR/input/$k\-big.in
	fi
done
