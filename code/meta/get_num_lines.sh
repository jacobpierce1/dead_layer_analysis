#!/bin/bash 
# find all non-peakdetect files with .py or .cpp suffix and count num lines.
cat $(find ../  -type f -not -path "*peakdetect*" -not -path "*#*" \
      -not -path "*lmfit*" -name "*.py" -o -name "*.cpp" -o -name "*.sh") | wc -l
