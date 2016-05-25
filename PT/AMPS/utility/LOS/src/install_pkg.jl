Pkg.add("DataFrames")
Pkg.add("HDF5")
Pkg.add("JLD")
#Pkg.update()

currentDir = pwd()
cd()
homeDir = pwd()
cd(currentDir)

fileName = ".juliarc.jl"
if !isfile(joinpath(homeDir, fileName))
  println(" - creating new .juliarc.jl file in home directory")
  fid = open(joinpath(homeDir, fileName), "w")
  write(fid, "push!(LOAD_PATH, pwd())\n")
  close(fid)
else
  println(" - checking .juliarc.jl for LOAD_PATH")
  hasPath = false
  fid = open(joinpath(homeDir, fileName), "r")
  while !eof(fid)
    line = readline(fid)
    if contains(line, "LOAD_PATH, pwd()")
      hasPath = true
    end
  end
  close(fid)
  if !hasPath
    println(" - appending pwd to LOAD_PAT")
    fid = open(joinpath(homeDir, fileName), "a")
    write(fid, "push!(LOAD_PATH, pwd())\n")
    close(fid)
  else
    info(".jularc.lj already has pwd in LOAD_PATH")
  end
end

