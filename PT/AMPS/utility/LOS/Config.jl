function compileLinux(tmpdir)
  previousDir = pwd()
  cd(joinpath(tmpdir, "lib"))
  run(`ar -x cspice.a`)
  run(`ar -x csupport.a`)
  objectFiles = AbstractString[]
  for fileName in readdir(".")
    if split(fileName, '.')[end] == "o"
      push!(objectFiles, fileName)
    end
  end
  run(`gcc -shared -fPIC -lm $objectFiles -o spice.so`)
  cd(previousDir)
  return objectFiles
end

function compileOSX(tmpdir)
  previousDir = pwd()
  cd(joinpath(tmpdir, "lib"))
  run(`ar -x cspice.a`)
  run(`ar -x csupport.a`)
  objectFiles = AbstractString[]
  for fileName in readdir(".")
    if split(fileName, '.')[end] == "o"
      push!(objectFiles, fileName)
    end
  end
  run(`gcc -dynamiclib -lm $objectFiles -o spice.dylib`)
  cd(previousDir)
  return objectFiles
end

function linux()
  global os
  os =  "linux"
end

function osx()
  global os
  os = "osx"
end


function updateSettings(keyword, newValue)
  mv(".userSettings.conf", ".userSettings.bkp")
  iFile = open(".userSettings.bkp", "r")
  oFile = open(".userSettings.conf", "w")
  foundKeyword = false
  while !eof(iFile)
    line = readline(iFile)
    if contains(line, keyword)
      if length(newValue) >= 1
        write(oFile, keyword, newValue, "\n")
      end
      foundKeyword = true
    else
      write(oFile, line)
    end
  end

  if (!foundKeyword) & (length(newValue) > 0)
    write(oFile, keyword, newValue, "\n")
  end

  close(iFile)
  close(oFile)
  rm(".userSettings.bkp")
end

function get_tmp_dir()
  keyword = "tmpDir:"
  iFile = open(".userSettings.conf", "r")
  while !eof(iFile)
    line = readline(iFile)
    if contains(line, keyword)
      value = string(bytestring(split(line, keyword)[2][1:end-1]))
      close(iFile)
      return value
    end
  end
  println(" - tmpDir not found in .userSettings.conf.")
  println(" - run julia Config.jl --tmpdir <path to tmpdir>")
  close(iFile)
  exit()
end


################################################################################
# start main
################################################################################

os = "operatingSystem"
@linux?  linux() : osx()

if lowercase(ARGS[1]) == "--tmpdir"
  rundir = ""
  try
    rundir = ARGS[2]
  catch
    rundir = pwd()
  end

  if isdir(rundir)
    println(" - tmpdir already exists")
    if (isdir(joinpath(rundir, "lib")) | isdir(joinpath(rundir, "input")))
      print(" - directory 'lib' or 'input' already exist, remove them? (y/n): ")
      answer = readline(STDIN)
      if contains(lowercase(answer), "y")
        try
          rm(joinpath(rundir, "lib"), recursive=true)
        catch
        end
        try
          rm(joinpath(rundir, "input"), recursive=true)
        catch
        end
        mkdir(joinpath(rundir, "lib"))
        mkdir(joinpath(rundir, "input"))
      end
    else
      println(" - create new directories 'lib' and 'input'")
      mkdir(joinpath(rundir, "lib"))
      mkdir(joinpath(rundir, "input"))
    end
  else
    mkdir(rundir)
    mkdir(joinpath(rundir, "lib"))
    mkdir(joinpath(rundir, "input"))
  end
  touch(".userSettings.conf")
  updateSettings("tmpDir:", rundir)

elseif lowercase(ARGS[1]) == "--spicelib"
  sharedLibPath = ARGS[2]
  tmpdir = get_tmp_dir()
  cp(joinpath(sharedLibPath, "cspice.a"),
     joinpath(tmpdir, "lib/cspice.a"),
     remove_destination=true)

  cp(joinpath(sharedLibPath, "csupport.a"),
     joinpath(tmpdir, "lib/csupport.a"),
     remove_destination=true)

  if os == "linux"
    objectFiles = compileLinux(tmpdir)
    updateSettings("spicelib:", joinpath(tmpdir, "lib/spice.so"))
  else
    objectFiles = compileOSX(tmpdir)
    updateSettings("spicelib:", joinpath(tmpdir, "lib/spice.dylib"))
  end
  previousDir = pwd()
  cd(joinpath(tmpdir, "lib"))
  run(`rm $objectFiles`)
  cd(previousDir)


elseif lowercase(ARGS[1]) == "--kernelfile"
  kernelFile = ARGS[2]
  updateSettings("kernelFile:", kernelFile)

elseif lowercase(ARGS[1]) == "--clib"
  clibFile = ARGS[2]
  tmpdir = get_tmp_dir()
  try
    cp(clibFile, joinpath(tmpdir, "lib/"*basename(clibFile)), remove_destination=true)
    updateSettings("clibFile:", joinpath(tmpdir, "lib/"*basename(clibFile)))
  catch
    println("Did not find c library! ")
    updateSettings("clibFile:", "")
  end

elseif lowercase(ARGS[1]) == "--noclib"
  updateSettings("clibFile:", "")

elseif lowercase(ARGS[1]) == "--meshfile"
  meshFile = ARGS[2]
  tmpdir = get_tmp_dir()
  cp(meshFile, joinpath(tmpdir, "input/"*basename(meshFile)), remove_destination=true)
  updateSettings("meshFile:", joinpath(tmpdir, "input/"*basename(meshFile)))

elseif lowercase(ARGS[1]) == "--meshfileshadow"
  meshFile = ARGS[2]
  tmpdir = get_tmp_dir()
  cp(meshFile, joinpath(tmpdir, "input/"*basename(meshFile)), remove_destination=true)
  cp(meshFile, joinpath(tmpdir, "input/"*basename(meshFile)), remove_destination=true)
  updateSettings("meshFileShadow:", joinpath(tmpdir, "input/"*basename(meshFile)))

elseif lowercase(ARGS[1]) == "--docheckshadow"
  doCheckShadow = ARGS[2]
  updateSettings("doCheckShadow:", doCheckShadow)

elseif lowercase(ARGS[1]) == "--datafile"
  dataFile = ARGS[2]
  tmpdir = get_tmp_dir()
  cp(dataFile, joinpath(tmpdir, "input/"*basename(dataFile)), remove_destination=true)
  updateSettings("dataFile:", joinpath(tmpdir, "input/"*basename(dataFile)))

elseif lowercase(ARGS[1]) == "--clean"
  tmpdir = get_tmp_dir()
  previousDir = pwd()
  cd(joinpath(tmpdir, "lib"))
  allFiles = readdir()
  if length(allFiles) > 0
    run(`rm $allFiles`)
  end
  cd(joinpath(tmpdir, "input"))
  allFiles = readdir()
  if length(allFiles) > 0
    run(`rm $allFiles`)
  end
  cd(previousDir)

elseif lowercase(ARGS[1]) == "--show"
  iFile = open(".userSettings.conf", "r")
  while !eof(iFile)
    print(" \u2764 " * readline(iFile))
  end
  close(iFile)
elseif lowercase(ARGS[1]) == "--help"
  print("say please: ")
  answ = readline(STDIN)
  if contains(answ, "please")
    println("Thanks! You have the following options to call Config.jl:")
    println("")
    println("--tmpdir          directory where files for runs are stored")
    println("--spicelib        directory to cspice.a and csupport.a")
    println("--kernelfile      spice kernel metafile")
    println("--clib            custom shared library to be used in LOS calculation")
    println("--noclib          do not use any shared c library")
    println("--meshfile        .ply file of the body surface mesh")
    println("--meshfileshadow  .ply file of the body surface mesh for shadow calc.")
    println("--docheckshadow   yes or no if shadow calculation is needed")
    println("--datafile        .full path to h5 AMPS output file")
    println("--clean           remove 'lib' and 'input' dirs in tmpfile")
    println("--help            show this message")
  end

end
