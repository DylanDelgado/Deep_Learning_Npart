#!/bin/bash

eng=3
dnum=800
jnum=8000

curdir=`pwd`

sed -i 's/RUNNAME=%NAME%/RUNNAME=u.'$eng'/' conf.mk
sed -i 's/%JNUM%/'$jnum'/'                  conf.mk
sed -i 's/@ENG@/'$eng'/'                    urqmd.cc
sed -i 's/@DNUM@/'$dnum'/'                  urqmd.cc
sed -i 's/@JNAME@/u.'$eng'/'                queue_daemon

cd $curdir
