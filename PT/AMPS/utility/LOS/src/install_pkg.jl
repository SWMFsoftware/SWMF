Pkg.add("DataFrames")
Pkg.add("HDF5")
Pkg.add("JLD")

if !isfile("~/.juliarc.jl")
  fid = open("~/.juliarc.jl", "w")
  write(fid, "push!(LOAD_PATH, pwd())\n")
  close(fid)
else
  hasPath = false
  fid = open("~/.juliarc.jl", "r")
  while !eof(fid)
    line = readline(fid)
    if contains(line, "LOAD_PATH, pwd()")
      hasPath = true
    end
  end
  close(fid)
  if !hasPath
    fid = open("~/.juliarc.jl", "a")
    write(fid, "push!(LOAD_PATH, pwd())\n")
    close(fid)
  end
end

    

  
