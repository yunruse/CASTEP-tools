#!/bin/bash -l

for f in cells/*.md; do
  printf "$(echo $f | cut -d/ -f2 | cut -d. -f1),"
  echo $(grep -e "T" $f | wc -l)
done
