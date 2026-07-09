#!/bin/bash

curdir=`pwd`
default_dist=./roots

if [ -z "$1" ]; then
	dist=$default_dist
else
	dist="$1"
fi

if [ -d $dist ]; then
	mv $dist $dist.`date -u '+%s'`
fi
mkdir $dist

count=1

for dir in Data/*_run*; do
	cd $dir
	echo $dir
	for rfile in *.root; do
		if ! [ -e "$rfile" ]; then
            echo $dir >> $curdir/noroot
			continue
		fi
		printf -v fname '%05d.root' $count
		((count=count+1))
		mv $rfile $curdir/$dist/$fname
	done
	cd $curdir
done
