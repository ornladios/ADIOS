#!/bin/bash

cd ~/gary
/home/liangz/local/mdtmstream/1.0.0-alpha.1/sbin/mdtm-ftp-server -data-interface 131.225.2.31 -password-file passfile -p 5501 -c server.conf -l mdtmftp.log -log-level all & export APP_PID=$!

cd /home/gary/demo/examples/C/global-array
python read_xgc.py -m -p 5565

kill -9 $APP_PID
