<?php  // If set, use passed value, otherwise use empty string
$parameter = array('margin', 'cmp', 'runname', 'plotstyle', 'plotfile', 'plottype');
foreach($parameter as $name) { $$name = isset($_GET[$name]) ? $_GET[$name] : ''; } ?>
