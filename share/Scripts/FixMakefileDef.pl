#!/usr/bin/perl -pi
s/^OS *=.*\n/"OS = ".`uname`/e;
s/^SWMF_ROOT *=.*\n/"SWMF_ROOT = ".`\/bin\/pwd`/e;
