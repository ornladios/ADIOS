#!/bin/bash

cd ~/gary
/home/liangz/local/mdtmstream/debug/mdtmftp/bin/mdtm-ftp-client -vb -p 4 file:///tmp/MdtmManPipes/ ftp://mdtmftp:123456@131.225.2.31:5501/tmp/ & export APP_PID=$!

cd /home/gary/demo/examples/C/global-array
./test_xgc

kill -9 $APP_PID
