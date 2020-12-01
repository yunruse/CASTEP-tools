#!/bin/bash -l

for f in cells/*.md; do
  perl -e "printf('%05.2f,', $(grep -e "T" $f | wc -l) / 2000)"
  echo $f | cut -d/ -f2 | cut -d. -f1
done
