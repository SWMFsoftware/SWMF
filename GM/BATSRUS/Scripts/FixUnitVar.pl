#!/usr/bin/perl -pi~
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# protect UnitUser_V
s/\bunituser_v\b/unituser\#v/gi;    

s/\bunituser_([a-z]+)\b/No2Io_V\(Unit\L\u$1\E_\)/gi;
s/\bunitsi_([a-z]+)\b/No2Si_V\(Unit\L\u$1\E_\)/gi;

# Use reverse transformation array instead of division
s/\/(\s*)No2Io_V/\*$1Io2No_V/g;
s/\/(\s*)No2Si_V/\*$1Si2No_V/g;

# Rename unitstr_tec and unistr_idl strings
s/\bunitstr_tec_([a-z]+)\b/NameTecUnit_V\(Unit\L\u$1\E_\)/gi;
s/\bunitstr_idl_([a-z]+)\b/NameIdlUnit_V\(Unit\L\u$1\E_\)/gi;

# recover UnitUser_V
s/unituser\#v/UnitUser_V/g;         
