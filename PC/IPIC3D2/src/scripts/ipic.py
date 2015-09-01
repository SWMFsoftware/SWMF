#!/usr/bin/env python

import os
import os.path
import sys
import subprocess
import socket # gethostname()
import re # regular expression
#from optparse import OptionParser
import getopt
# http://docs.python.org/2/library/collections.html#collections.deque
from collections import deque # double-ended queue
import inspect
#
# useful documentation:
#
# http://effbot.org/zone/python-list.htm
# http://pymotw.com/2/subprocess/
# http://stackoverflow.com/questions/3777301/how-to-call-a-shell-script-from-python-code

def getdims(inputfile):
    # extract dimensions from intput file
    dims = [1, 1, 1]
    pattern = re.compile(r'^\s*([\w]+)\s*=\s*([\w]+)')
    f = open(inputfile)
    for line in f:
        # key, value = line.split('=')
        #pattern.findall(line)
        match = re.search(pattern, line)
        if match:
            var = match.group(1)
            val = match.group(2)
            if var == 'XLEN':
                dims[0]=int(val)
            elif var == 'YLEN':
                dims[1] = int(val)
            elif var == 'ZLEN':
                dims[2] = int(val)
    f.close()
    return dims

def lineno():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno

def capture_shell_command_output(command):
  return os.popen(command).read()

def issue_command(command):
  if(show):
    print ' '.join(command)
  else:
    print '+', ' '.join(command)
    subprocess.call(command);

def issue_shell_command(command):
  if(show):
    print command
  else:
    print '+', command
    os.system(command)

def construct_run_command(args,mpirun):

    # convert from deque to list for getopts
    args = list(args)

    # set default values
    num_nodes = 1
    num_max_threads = 1
    num_threads_per_node = 1
    outputdir = 'data'
    inputfile = 'parameters.inp'
    hostname = ''
    # these parameters should be specified in ipic module
    ipic_env = os.getenv('IPIC_ENV','default');
    ipic_runenv = os.getenv('IPIC_RUNENV',ipic_env);
    if ipic_runenv == 'miclogin.xeon' or ipic_runenv == 'miclogin.mic':
      if ipic_runenv == 'miclogin.xeon':
        mpirun = 'mpiexec.hydra' # is this line needed?
        # calculate number of threads per process
        # + could extract this stuff from /proc/cpuinfo
        #   via commands such as:
        #   grep processor /proc/cpuinfo  | wc # gives 32
        #   grep -m 1 siblings /proc/cpuinfo # gives 16
        #   grep -m 1 cores /proc/cpuinfo # gives 8
        # + key:
        #   - siblings reside on same processor (16=8*2)
        #   - cores is number of cores per processor
        #   - "processor" = hardware thread
        num_nodes = 1
        num_processors_per_node = 2
        num_cores_per_processor = 8
        num_threads_per_core = 2
        num_threads_per_node = (
          num_threads_per_core *
          num_cores_per_processor *
          num_processors_per_node)
      elif ipic_runenv == 'miclogin.mic':
        mpirun = 'mpiexec.hydra'
        # calculate number of threads per process
        # - could use ssh to extract this stuff from /proc/cpuinfo
        #   on the machine we will run on via e.g.
        #   ssh $hostname grep -m 1 cores /proc/cpuinfo | awk '{print $4}'61
        num_nodes = 1
        num_processors_per_node = 1
        hostname = socket.gethostname()
        num_cores_per_processor = 60
        if hostname=='knc1':
          num_cores_per_processor = 61
          # 0,1,2,3 all exist
          micnum = 3
        elif hostname=='knc2':
          num_cores_per_processor = 57
          micnum = 0
        else:
          print 'hostname', hostname, 'not supported'
          sys.exit(-1)
        num_threads_per_core = 4
        num_threads_per_node = (
          num_threads_per_core *
          num_cores_per_processor *
          num_processors_per_node)
        #
        hostname = hostname + '-mic' + str(micnum)

    # now that the default mpirun has been determined,
    # allow $IPIC_MPIRUN to override it
    mpirun = os.getenv('IPIC_MPIRUN',mpirun);

    num_threads_is_given_by_user = 0
    try:
      opts, args = getopt.getopt(args, 'i:o:t:h:', \
        ['input=', 'outputdir=', 'threads=', 'host='])
    except getopt.GetoptError, e:
      if e.opt == 'h' and 'requires argument' in e.msg:
        print 'ERROR: -h requires input filename'
      elif e.opt == 'i' and 'requires argument' in e.msg:
        print 'ERROR: -i requires input filename'
      elif e.opt == 'o' and 'requires argument' in e.msg:
        print 'ERROR: -o requires directory name'
      elif e.opt == 't' and 'requires argument' in e.msg:
        print 'ERROR: -t requires max number of threads'
      else:
        usage_error()

    for o, a in opts:
        if o in ("-h", "--host"):
          hostname = a
        elif o in ("-i", "--input"):
          inputfile = a
        elif o in ("-o", "--outputdir"):
          outputdir = a
          print 'ERROR: -o is not yet supported'
          sys.exit(1)
        elif o in ("-t", "--threads"):
          num_threads_is_given_by_user = 1
          num_max_threads = int(a)
        #else:
        #  assert False, "unhandled option"

    if len(args)!=0:
      usage_error()

    # determine num_procs
    dims = getdims(inputfile)
    XLEN = dims[0]
    YLEN = dims[1]
    ZLEN = dims[2]
    num_procs = XLEN*YLEN*ZLEN
    num_procs_per_node = num_procs/num_nodes
    num_threads_per_proc = num_threads_per_node/num_procs_per_node

    if not num_threads_is_given_by_user:
      # rounding down is the correct behavior
      num_max_threads = int(num_threads_per_proc)

    arguments = ['./iPic3D', outputdir+"/parameters.inp"];
    options = ['-np', str(num_procs)]
    if hostname!="":
        options.extend(['-host', hostname])

    if num_max_threads > 1:
        omp_string = 'OMP_NUM_THREADS=' + str(num_max_threads)
        omp = ['-env', omp_string]
        options.extend(omp)

    command = [mpirun]
    command.extend(options)
    command.extend(arguments)
    # add any additional arguments given by user
    #command.extend(args)
    return outputdir, inputfile, command

def ipic_run(args,mpirun):
    # clear the data directory
    # (currently by H5hut, i.e. when setting WriteMethod=Parallel)
    outputdir, inputfile, command = construct_run_command(args,mpirun)
    # issue_command(['rm', '-rf', outputdir])
    issue_shell_command('rm -rf ' + outputdir + '/*')
    issue_command(['mkdir', '-p', outputdir])
    issue_command(['cp', inputfile, outputdir+"/parameters.inp"])
    issue_command(command)

def ipic_make_data():
    # create data subdirectory
    create_data_command = '''mkdir -p data''';
    issue_shell_command(create_data_command)

def get_inputfile(args):
    # convert from deque to list for getopts
    args = list(args)
    # set default inputfile
    inputfile = 'src/inputfiles/GEM.inp'
    try:
      opts, args = getopt.getopt(args, 'i:', ['input='])
    except getopt.GetoptError, e:
      if e.opt == 'i' and 'requires argument' in e.msg:
        print 'ERROR: -i requires input filename'
      else:
        usage_error()
    #
    for o, a in opts:
        if o in ("-i", "--input"):
          inputfile = a
        #else:
        #  assert False, "unhandled option"
    #
    if len(args)>=2:
      usage_error()
    return inputfile, args

def ipic_cmake(args):

    # get the input file
    inputfile, args = get_inputfile(args)

    # extract the values of any options from args
    #system, parhdf5, args = get_cmake_options(args)

    # make src a link to the code
    numargs = len(args)
    if numargs==0:
      sourcedir = os.getenv('IPIC_HOME','..');
    else:
      args = deque(args)
      sourcedir = deque.popleft(args)

    if sourcedir!='src':
      rm_command = ['rm', '-f', 'src'];
      issue_command(rm_command);
      ln_command = ['ln', '-s', str(sourcedir), 'src'];
      issue_command(ln_command)

    ipic_env = os.getenv('IPIC_ENV','')
    if ipic_env != '':
      ipic_home = os.getenv('IPIC_HOME','')
      if ipic_home != '':
        jobscript_file = ipic_home+'/env/'+ipic_env+'/job.sh'
        if os.path.isfile(jobscript_file):
          command = ['cp', jobscript_file, '.'];
          issue_command(command)
      
    # copy input file to current directory
    copy_input_file_command = ['cp', inputfile, './parameters.inp']
    issue_command(copy_input_file_command)

    ipic_make_data();

    # invoke cmake 
    cmake = os.getenv('IPIC_CMAKE','cmake');
    cmake_command = [cmake];

    # add arguments particular to intel
    #cmake_command.extend(['-DCMAKE_CXX_FLAGS=-DMPICH_IGNORE_CXX_SEEK'])

    # add rest of args to end of cmake command

    # issue the command
    cmake_command.extend(['src'])
    cmake_command.extend(args)
    # a string is a list of individual characters,
    # so we append rather than extend the list
    cmake_command.append(os.getenv('IPIC_CMAKE_ARGS',''))
    # by spawning a subshell we allow the content of
    # IPIC_CMAKE_ARGS to be interpreted by shell
    issue_shell_command(' '.join(cmake_command))

def ipic_make(args):

    # invoke make 
    make_command = ['time', 'make', 'VERBOSE=1'];
    make_command.extend(args)
    #issue_shell_command(' '.join(make_command))
    issue_command(make_command)

def ipic_ctags(args):
    # create tags file using ctags
    create_tags_command = \
        '''find . -name '*.cpp' -or -name '*.h' | grep -v unused | xargs ctags --extra=+qf'''
    issue_shell_command(create_tags_command)
    # sort tags file
    sort_tags_command = '''LC_ALL=C sort -u tags -o tags'''
    issue_shell_command(sort_tags_command)

def ipic_eval_shell(command, args):
    ipic_command = ['ipic-'+command]
    ipic_command.extend(args)
    issue_command(ipic_command)

#def ipic_eval_shell(command, args):
#    showcommand = ['ipic-show-'+command]
#    showcommand.extend(args)
#    shell_command = " ".join(showcommand)
#    generated_command = capture_shell_command_output(shell_command)
#    # generated_command = capture_command_output(showcommand)
#    issue_shell_command(generated_command)

def ipic_basic_help():
    print '''
  To build, you can use:
  
    mkdir build
    cd build
   ''', progname, '''cmake /path/to/ipic3d
    make # or "make -j" to compile in parallel
  
  Then to run the code, use:
  
    ipic run

  If you prefer, use e.g. "ipic show run" to see the shell commands
  that will be executed and then execute them directly yourself.
  
  Minor subcommands:

    ''', progname, '''help ctags   # create ctags file to navigate code
    ''', progname, '''help mic     # help for running on mic
    ''', progname, '''help deep    # help for running on deep

  Major subcommands:

    ''', progname, '''help show    # show what a command would do
    ''', progname, '''help cmake   # execute cmake and create subdirectories
    ''', progname, '''help run     # execute iPic3D
  '''

def ipic_help_show(args):
    print '''
  ''', progname, '''show [command]

    show the shell command that would be executed by
      ipic [command]
    '''

def ipic_help_run(args):
    print '''
 ''', progname, '''run [options]

    run iPic3D with appropriate arguments.

    options:
    -t <num_max_threads>: set maximum number of threads per process
       (default is 1 unless -s <mic|xeon> is set)
    -i <inputfile>: set input file (default is "src/inputfiles/GEM.inp")
    -o <outputdir>: set output directory (default is "data")
    -h <host>: spawn processes on specified host
    '''

def ipic_help_mic(args):
    print '''
  Commands to run the code on miclogin:

    ssh miclogin
    # check out the code
    git clone https://github.com/alecjohnson/iPic3D.git ipic3d

    # put the following lines in your .bashrc:
    #
    export IPIC_HOME=$HOME/ipic3d
    alias ipic="$IPIC_HOME/scripts/ipic"                                          
    module use $IPIC_HOME/sysenv/miclogin                                         

    # log into the (confusingly named) Xeon host processor:
    ssh knc1 # or ssh knc2

    # compile for and run on xeon
    #
    module purge
    # use parallel HDF5 for I/O
    module load ipic-parallel
    # use non-parallel HDF5 for I/O
    #module load ipic
    rm -rf build; mkdir build; cd build
    ipic cmake $IPIC_HOME
    make -j VERBOSE=1
    ipic run

    # compile for and run on MIC
    #
    module purge
    # use parallel HDF5 for I/O
    module load ipic-mic-parallel
    # use non-parallel HDF5 for I/O
    #module load ipic-mic
    rm -rf build.mic; mkdir build.mic; cd build.mic
    ipic cmake $IPIC_HOME
    make -j VERBOSE=1

  To show what a command will do, use e.g.:

    ipic show [command]
  
  See also:
    ''', progname, '''help
    ''', progname, '''help deep
    '''

def ipic_help_deep(args):
    print '''
  # put the following lines in your ~/.bashrc:
  #
  export IPIC_HOME=$HOME/ipic3d
  module use $IPIC_HOME/sysenv/deep

  # load the ipic module
  module load ipic

  # compile and run
  #
  mkdir build; cd build
  ipic cmake
  ipic make
  ipic run

  # If you want to use parallel I/O,
  # look at the notes in $IPIC_HOME/README
  # and $IPIC_HOME/env/deep/ipic-parallel
    '''

def ipic_help_cmake(args):
    print '''
  ''', progname, '''cmake [sourcedir]

     [sourcedir]: the source code directory; by default ${IPIC_HOME:-..}
  '''

def ipic_help_ctags(args):
    print '''
  Make sure that you are in the source code directory
  and then run

    ''', progname, '''ctags

  This will generate, display, and execute a ctags command
  that creates a file named "tags", usable by the vim editor:

    vim -t main # open editor at definition of main
    vim -c "help tag" # get help on using tag files.
    '''

def ipic_help_git(args):
    print '''
    ### This stub gives examples of git commands ###

    # show branch information
    git branch -avv
    # examining the .git directory reveals a wealth of information, e.g.:
    cat .git/config
    # with --stat all files checked in are displayed.
    git log --stat
    # for the following I just do "git tree" (see .gitconfig below):
    git log --oneline --decorate --graph --branches --source
    git status # shows file statuses
    git remote -v # show remote repositories
    # show commits in chronological order.
    git reflog
    # git reflog is useful to get the sha-1 hash of a commit
    # that you recently made and whose branch you accidentally
    # deleted, making it no longer reachable.  Note that
    # each snapshot that you commit should stay in its local
    # repository for 90 days before being garbage collected
    # unless you do something like "git gc".  See also
    # http://gitready.com/advanced/2009/01/17/restoring-lost-commits.html
    #
    # show file
    git show mybranch:myfile
    eg cat myfile # slightly nicer than git show
    # show who checked in what line when under what commit.
    git blame myfile # on current branch
    git blame amaya-library iPic3D.cpp

  for modification:

    # initialize a repository
    mkdir localrepository; cd localrepository; git init
    # creating/removing remote:
    git remote add myremote  https://github.com/alecjohnson/iPic3D.git
    git remote rm myremote  
    # get all branches and their filesystem snapshots
    # from myremote that are not already in localrepository
    git fetch myremote 
    # check in mods
    git stage myfile
    git rm oldfile
    git commit
    # modify a commit message
    git commit --amend
    # create a branch and check it out
    git checkout -b newbranch
    # push branch to server
    eg push --branch newbranch myremote
    # pull changes from server into current branch
    git pull
    # delete a branch on server (!):
    git push myremote --delete mybranch

  # example of global configuration file:

    $ cat ~/.gitconfig
    [user]
    name = eajohnson
    email = e.alec.johnson@gmail.com
    [alias]
    tree = log --oneline --decorate --graph --branches --source
    undo-commit = reset --soft HEAD~1
    '''

def ipic_help(args):
    if len(args) == 0:
      ipic_basic_help()
      sys.exit()
    
    command = deque.popleft(args)
    if command == "show":
      ipic_help_show(args)
    elif command == "run":
      ipic_help_run(args)
    elif command == "mic":
      ipic_help_mic(args)
    elif command == "deep":
      ipic_help_deep(args)
    elif command == "cmake":
      ipic_help_cmake(args)
    elif command == "ctags":
      ipic_help_ctags(args)
    elif command == "git":
      ipic_help_git(args)
    else:
        print "ipic help", command, "is not supported"
        sys.exit(-1)

def lprint(message):
    theline = inspect.currentframe().f_back.f_lineno
    print ' line', str(theline), ': ', message

def invalid_value_error(message, value):
    theline = inspect.currentframe().f_back.f_lineno
    print 'ERROR, line', str(theline), ':', message, value
    sys.exit(-1)

def error(message):
    theline = inspect.currentframe().f_back.f_lineno
    print 'ERROR, line', str(theline), ':', message
    sys.exit(-1)

def usage_message():
    print '''
  usage: ''', progname, ''' [show] <command>

  Available commands:
    ''', progname, '''help
    ''', progname, '''show
    ''', progname, '''cmake
    ''', progname, '''run
      '''

def usage_error():
  theline = inspect.currentframe().f_back.f_lineno
  print '  usage_error() called from ipic.py line ', str(theline)
  usage_message()
  sys.exit(-1)

def ipic_command(argv1):

    # it might be better to use the argparse module rather than getopt,
    # but unfortunately argparse is only available beginning with python 2.7
    # and most HPC platforms seem to have python 2.6 installed.
    # optparse has been deprecated and does not seem to be in python 3;
    # note, however, that argparse was initially an extension of optparse
    # before giving up on backward compatibility.
    #
    try:
      opts, args = getopt.getopt(argv1, 'h', ['help'])
    except getopt.GetoptError, e:
      usage_error()

    for o, a in opts:
        if o in ("-h", "--help"):
          usage_error()
        #else:
        #  assert False, "unhandled option"

    #print "type(args)", type(args)

    numargs = len(args)
    if numargs==0:
      usage_error()

    args = deque(args)
    command = deque.popleft(args)

    if command == "help":
        ipic_help(args)
    elif command == "ctags":
        ipic_ctags(args)
    elif command == "cmake":
        ipic_cmake(args)
    elif command == "make":
        ipic_make(args)
    elif command == "run":
        ipic_run(args, 'mpirun')
    elif command == "exec":
        ipic_run(args, 'mpiexec')
    elif command == "setstripe" \
      or command == "getstripe":
        ipic_eval_shell(command, args)
    else:
        print progname, command, "is not supported. Try: ipic help"
        sys.exit(-1)

    #print os.path.basename(__file__)
    #print os.path.dirname(__file__)

def main():

    global progname
    progname = os.path.basename(sys.argv[0])
    global dirname
    dirname = os.path.dirname(sys.argv[0])
    global show
    show=0

    argv1 = sys.argv[1:]
    if len(argv1)==0:
      usage_message()
      sys.exit()

    if argv1[0]=='show':
      show=1
      argv1=argv1[1:]

    ipic_command(argv1)

if __name__ == '__main__':
    main()

