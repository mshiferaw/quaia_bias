#!/bin/bash

for j in `seq 5749150 5749365` ; do
    scancel -u mahlet $j
    echo $j
done