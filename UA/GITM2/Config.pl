#!/usr/bin/perl -s

if ($Earth) {$earth=1;}
if ($EARTH) {$earth=1;}

if ($Mars) {$mars=1;}
if ($MARS) {$mars=1;}

if (!$mars && !$earth) {
    print "Must specify either -earth or -mars\n";
    $h = 1;
}

if ($h) { print "Help!!\n"; }


if ($earth) { 

    print "Configuring GITM for Earth!!\n"; 

    $command = "rm -f src/ModPlanet.f90 src/ModPlanet.o";
    if ($v) { print " -> $command,\n";}
    system $command;
    $command = "cd src ; ln -s ModEarth.f90 ModPlanet.f90";
    if ($v) { print " -> $command,\n";}
    system $command;

    $command = "rm -f src/planet.f90";
    if ($v) { print " -> $command,\n";}
    system $command;
    $command = "cd src ; ln -s Earth.f90 planet.f90";
    if ($v) { print " -> $command,\n";}
    system $command;

    foreach $file (glob("src/*.Earth.f90")) {
	if ($file =~ /src\/(.*).Earth.f90/) {
	    if ($v) { print "linking file -> $1\n";}
	    $command = "rm -f src/$1.f90";
	    if ($v) { print " -> $command,\n";}
	    system $command;
	    $command = "cd src ; ln -s $1.Earth.f90 $1.f90";
	    if ($v) { print " -> $command,\n";}
	    system $command;
	}
    }

}

if ($mars) { 

    print "Configuring GITM for Mars!!\n"; 

    $command = "rm -f src/ModPlanet.f90 src/ModPlanet.o";
    if ($v) { print " -> $command,\n";}
    system $command;
    $command = "cd src ; ln -s ModMars.f90 ModPlanet.f90";
    if ($v) { print " -> $command,\n";}
    system $command;

    $command = "rm -f src/planet.f90";
    if ($v) { print " -> $command,\n";}
    system $command;
    $command = "cd src ; ln -s Mars.f90 planet.f90";
    if ($v) { print " -> $command,\n";}
    system $command;

    foreach $file (glob("src/*.Mars.f90")) {
	if ($file =~ /src\/(.*).Mars.f90/) {
	    if ($v) { print "linking file -> $1\n";}
	    $command = "rm -f src/$1.f90";
	    if ($v) { print " -> $command,\n";}
	    system $command;
	    $command = "cd src ; ln -s $1.Mars.f90 $1.f90";
	    if ($v) { print " -> $command,\n";}
	    system $command;
	}
    }

}

