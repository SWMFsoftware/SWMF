#!/usr/bin/perl -p
# Usage:   cat_sat_log.pl logfile1 logfile2 .... > log_all.log
# Example: cat_sat_log.pl log_*.log > log_all.log
$_ = '' unless ($. == 1 or $. == 2 or /^[\s\d\.eEdD\+\-]+$/);
