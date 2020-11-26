#!/bin/bash

while [[ -n $(qstat) ]]; do
  clear
  qwait
  sleep 20
done
echo -e "\a"
clear
echo "Queue is clear"
