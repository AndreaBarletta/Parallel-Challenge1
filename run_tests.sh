#!/bin/bash

arg=10000

for ((i=1; i<=15; i++)); do
  ./mergesort-ser $arg
  ./mergesort-par $arg
  arg=$((arg * 2))
done