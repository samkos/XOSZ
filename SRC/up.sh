#!/bin/sh
perfdmf_configure --create-default

# Just hit enter
perfexplorer_configure  < in
echo "Perfexplorer configured..."

for i in  8 16 32 64 128 256; do
  echo "Uploading data for lm=nm=$i"
  perfdmf_loadtrial -a ZEPHYR -x Poincare -n $i "../RES/$i/scorep/profile.cubex"
done;