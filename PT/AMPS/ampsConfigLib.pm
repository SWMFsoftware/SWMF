#$Id$

use strict;
use warnings;

use POSIX qw(strftime);
use List::Util qw(first);

package ampsConfigLib;
our @EXPORT=qw($WorkingSourceDirectory);

my $WorkingSourceDirectory;


sub printwd {
  print "$ampsConfigLib::WorkingSourceDirectory\n";
}

#=============================== Change a value of a variable in the code  =============================
sub ChangeValueOfVariable {
  my $Variable=$_[0];
  my $Value=$_[1];
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
  
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File");
  
  
  foreach (@FileContent) {
    if ($_=~/($Variable)/) {
#      $Variable=~s/\\//g;
#      $_=$Variable."=".$Value.";\n";

      my $t=$Variable;
      
      $t=~s/\\//g;
      $_=$t."=".$Value.";\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT);
}

#=============================== Change a definition of a macro in a source code  =============================
sub RedefineMacro {
  my $Macro=$_[0];
  my $Value=$_[1];
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
  
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";
  
  foreach (@FileContent) {
    if ($_=~/\#define $Macro /) {
#      print "\#define $Macro $Value\n";
      $_="\#define $Macro $Value\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT);
}

sub RedefineMacroFunction {
  my $Macro=$_[0];
  my $Value=$_[1];
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
  
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";
  
  foreach (@FileContent) {
    if ($_=~/\#define $Macro/) {
#      print "\#define $Macro$Value\n";
      $_="\#define $Macro$Value\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT);
}

sub AddMacro {
  my $Macro=$_[0];
  my $Value=$_[1];
  my $File=$_[2];
  
  open (FILEIN,">>$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";   

  print FILEIN  "//The following changes are added while reading input file. [".POSIX::strftime("%m/%d/%Y", localtime)."]\n";
  print FILEIN "#undef $Macro\n#define $Macro $Value\n\n";
  close (FILEIN);
}
  
#=============================== Change a value of a array in the code  =============================
sub ChangeValueOfArray {
  my $Variable=$_[0];
  my @Value=@{$_[1]};
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
 
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File");
  
  foreach (@FileContent) {
    if ($_=~m/($Variable)/) {
      my $t=$Variable;
      
      $t=~s/\\//g;
      $_=$t."={";
      
      for (my $i=0;$i<1+$#Value;$i++) {
        $_=$_."$Value[$i]";
        
        if ($i ne $#Value) {
          $_=$_.",";
        }
      }
      
      $_=$_."};\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT); 
}

#===============================  Get the first appearence of a string $_[0] in the array $_[1] ==
sub GetElementNumber {
  my $str=$_[0];
  my @List=@{$_[1]};
  
  my $nspec=0;
  my $found=0;
  
  foreach (@List) {
    if ($List[$nspec] eq $str) {
      $found=1; #the species has been found in the list of those used in the simulation
      last;
    }
    else {
      $nspec++;
    }
  }
  
  if ($found == 0) {
    $nspec=-1;
  }

  return $nspec;
}

#===============================  Substitute a line in a source file =============================
sub SubstituteCodeLine {
  my $oldLinstKey=$_[0];
  my $newLine=$_[1];
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
 
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File");
  
  foreach (@FileContent) {
    if ($_=~m/($oldLinstKey)/) {
      $_=$newLine."\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT);
}

#===============================  Insert a line into a file after a marker ==========================
sub AddLineAfterMarker2File {
  my $NewLine=$_[0];
  my $Marker=$_[1];
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
 
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File");
  
  foreach (@FileContent) {
    print FILEOUT "$_";
    
    if ($_=~m/($Marker)/) {
      print FILEOUT  "//The following changes are added while reading input file. [".POSIX::strftime("%m/%d/%Y", localtime)."]\n";
      print FILEOUT "$NewLine\n";
    }
  }
  
  close (FILEOUT);
}  


#===============================  RECURSIVE SUBSTITUDE OF A STRING IN THE SOURCE DIRECTORY ==========
sub RecursiveSubstitute {
  my $init=$_[0];
  my $final=$_[1];
  my $dir=$_[2];
  
  my $fname;
  my @FileList;
  
  opendir(DIR,$dir) or die "Cannot open directory $dir\n";
  @FileList=readdir(DIR);
  closedir(DIR);
  
  foreach $fname (@FileList) {
    #check if '$fname' is a directory
    
    if (($fname eq '.')||($fname eq '..')||($fname eq 'CVS')||($fname eq 'RCS')) {
      next;
    }
    
    if (-d "$dir/$fname") {
      #fname is a directory -> call the function recursevely
      RecursiveSubstitute($init,$final,"$dir/$fname");
    }
    else {
      #fname is a file -> make the substitution
      my @lines;
      
      open (FILE,"<$dir/$fname") || die "Cannot open file $dir/$fname\n";
      @lines=<FILE>;
      close (FILE);
      
      foreach (@lines) {
        $_=~s/$init/$final/g;
      }
      
      open (FILE,">$dir/$fname") || die "Cannot open file $dir/$fname\n";
      print FILE @lines;
      close (FILE);
    }
  }
}

sub AddLine2File {
  my $fname=$_[1];
  my $newline=$_[0];
  
  open (EDITFILE,">>$ampsConfigLib::WorkingSourceDirectory/$fname") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$fname\n";

  print EDITFILE  "//The following changes are added while reading input file. [".POSIX::strftime("%m/%d/%Y", localtime)."]\n";
  print EDITFILE "$newline\n";
  close(EDITFILE);
}

1;
