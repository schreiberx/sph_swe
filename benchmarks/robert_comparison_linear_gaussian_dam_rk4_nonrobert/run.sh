#! /bin/bash

BASEDIR="`pwd`"
rm -f "$BASEDIR/prog_*"

cd ../../

make clean
make release


# dimensional
#./build/sh_example T64 P0 M2 R$((2*60)) E$((5*60*60)) V0.0

# h0=g=f=1

# Stable time step size for T64
TS=0.001

# unstable for explicit RK4 time stepping
#TS=0.01

OTS=0.01

./build/sh_example T64 P0 M6 R$TS O$OTS E5 V0.0 F0 S0


mv prog_* "$BASEDIR"

cd "$BASEDIR"
../../plot.py prog_h_*.csv
../../create_mp4.sh prog_h out_prog_h.mp4
