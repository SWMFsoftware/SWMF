using HDF5, JLD
using PyPlot
using Spice
using Instrument

include("io.jl")
include("octree.jl")
include("raytrace.jl")

function interpolate(nVars, coords, oct)
  nPoints = size(coords, 2)
  data = zeros(Float64, nVars)
  myPoint = zeros(Float64, 3)
  result = zeros(Float64, nVars, nPoints)
  norms = Float64[]
  for i=1:nPoints
    for k=1:3
      myPoint[k] = coords[k,i]
    end

    didFindCell, myCell = cell_containing_point(oct, myPoint)
    if (didFindCell == true)
      triLinearInterpolation!(myCell, myPoint, data, 0.0, 0.0, 0.0)
    end

    for j=1:nVars
      result[j,i] = data[j]
      data[j] = 0.0
    end
  end
  return result
end

function save_results(nVars, result, coords, fileName="interp_output.dat")
  nPoints = size(result, 2)
  oFile = open(fileName, "w")
  for i=1:nPoints
    for k=1:3
      @printf(oFile, "%.9e ", coords[k,i])
    end
    for k=1:nVars
      @printf(oFile, "%.5e ", result[k,i])
    end
    @printf(oFile, "\n")
  end
end

################################################################################
# START main
################################################################################

global const clib = parseUserFile("clibFile:")
@show(clib)

coordFileName = ARGS[1]
metaFile = parseUserFile("kernelFile:")
try
  furnsh(metaFile)
catch e
  println("spice error")
  println(e)
  exit()
end
const fileName = parseUserFile("dataFile:")
const filePath = dirname(fileName)
fileNameExtension = split(basename(fileName), ".")[end]
fileNameBase = basename(fileName)[1:end-(length(fileNameExtension)+1)]


oct, nVars = build_octree(filePath, fileNameBase)
dummyCell = Cell(zeros(3),
                 zeros(3),
                 zeros(8,3),
                 0.0,
                 zeros(2,2),
                 false,
                 Triangle[],
                 0)
oct.cells[1] = dummyCell

println(" - nVars : ", nVars)
coords = load_user_coordinates(coordFileName, 6)
coords = coords[4:6,:]

@time result = interpolate(nVars, coords, oct)
@time save_results(nVars, result, coords, "interp_output.H2O.dat")
