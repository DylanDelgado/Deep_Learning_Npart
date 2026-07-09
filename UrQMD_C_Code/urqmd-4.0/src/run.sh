#!/bin/bash
#$ -cwd
#$ -o job.log 
#$ -e job.error
#$ -j y

echo 1 > tmpl.pid
./convert
echo $status > tmpl.pid
