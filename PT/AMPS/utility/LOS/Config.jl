function compileLinux()
  cd("lib")
  run(`ar -x cspice.a`)
  run(`ar -x csupport.a`)
  objectFiles = AbstractString[]
  for fileName in readdir(".")
    if split(fileName, '.')[end] == "o"
      push!(objectFiles, fileName)
    end
  end
  run(`gcc -shared -fPIC -lm $objectFiles -o spice.so`)
  cd("..")
  return objectFiles
end

function compileOSX()
  run(`ar -x cspice.a`)
  run(`ar -x csupport.a`)
  objectFiles = AbstractString[]
  for fileName in readdir(".")
    if split(fileName, '.')[end] == "o"
      push!(objectFiles, fileName)
    end
  end
  run(`gcc -dynamiclib -lm $objectFiles -o spice.dylib`)
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
  cd("input")
  mv("userSettings.conf", "userSettings.bkp")
  iFile = open("userSettings.bkp", "r")
  oFile = open("userSettings.conf", "w")
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
  rm("userSettings.bkp")
  cd("..")
end


################################################################################
# start main
################################################################################

os = "operatingSystem"
@linux?  linux() : osx()

if lowercase(ARGS[1]) == "--spicelib"
  sharedLibPath = ARGS[2]
  cp(joinpath(sharedLibPath, "cspice.a"), "lib/cspice.a", remove_destination=true)
  cp(joinpath(sharedLibPath, "csupport.a"), "lib/csupport.a", remove_destination=true)
  if os == "linux"
    objectFiles = compileLinux()
    updateSettings("spicelib:", "../lib/spice.so")
  else
    objectFiles = compileOSX()
    updateSettings("spicelib:", "../lib/spice.dylib")
  end
  cd("lib")
  run(`rm $objectFiles`)
  cd("..")

elseif lowercase(ARGS[1]) == "--kernelfile"
  kernelFile = ARGS[2]
  updateSettings("kernelFile:", kernelFile)

elseif lowercase(ARGS[1]) == "--clib"
  clibFile = ARGS[2]
  try
    cp(clibFile, "lib/"*basename(clibFile), remove_destination=true)
    updateSettings("clibFile:", "lib/"*basename(clibFile))
  catch
    println("Did not find c library! ")
    updateSettings("clibFile:", "")
  end

elseif lowercase(ARGS[1]) == "--noclib"
  updateSettings("clibFile:", "")

elseif lowercase(ARGS[1]) == "--meshfile"
  meshFile = ARGS[2]
  cp(meshFile, "input/"*basename(meshFile), remove_destination=true)
  updateSettings("meshFile:", "../input/"*basename(meshFile))

elseif lowercase(ARGS[1]) == "--meshfileshadow"
  meshFile = ARGS[2]
  cp(meshFile, "input/"*basename(meshFile), remove_destination=true)
  updateSettings("meshFileShadow:", "../input/"*basename(meshFile))

elseif lowercase(ARGS[1]) == "--docheckshadow"
  doCheckShadow = ARGS[2]
  updateSettings("doCheckShadow:", doCheckShadow)

elseif lowercase(ARGS[1]) == "--datafile"
  dataFile = ARGS[2]
  cp(dataFile, "input/"*basename(dataFile), remove_destination=true)
  updateSettings("dataFile:", "../input/"*basename(dataFile))

elseif lowercase(ARGS[1]) == "--clean"
  cd("lib")
  allFiles = readdir()
  if length(allFiles) > 0
    run(`rm $allFiles`)
  end
  cd("../input")
  allFiles = readdir()
  if length(allFiles) > 0
    run(`rm $allFiles`)
  end

end
