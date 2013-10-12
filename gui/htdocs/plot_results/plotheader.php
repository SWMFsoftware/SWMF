<?php 
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf 
?>
<?php  // If set, use passed value, otherwise use empty string
$parameter = array('margin', 'cmp', 'runname', 'plotstyle', 'plotfile', 'plottype');
foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>
