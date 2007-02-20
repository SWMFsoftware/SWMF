#!/usr/bin/perl -pi~xs
s/\bunituser_v\b/unituser\#v/gi;    # protect UnitUser_V
s/\bunituser_([a-z]+)\b/No2Io_V(Unit\L\u$1\E_\)/gi;
s/\bunitsi_([a-z]+)\b/No2Si_V(Unit\L\u$1\E_\)/gi;
s/unituser\#v/UnitUser_V/g;         # recover UnitUser_V
