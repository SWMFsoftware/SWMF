using HDF5, JLD
include("io.jl")
include("octree.jl")
# this script has to be applied for all new 3D cases with amps.
# it reads the original amps output file and saves some content
# of it into binary hdf5 files for quicker loading when doing
# the actual line of sight calculations with newLOS2.jl

# This script has to be run only once and needs the filename
# of the AMPS output as only input parameter.
# example run command:
# julia prepareAmpsData.jl /www//www/ices/Data/Coma/DSMC/CG_3.5_au_50/CG_3.5_au_50.H2O.dat

fileNameExtension = split(basename(ARGS[1]), ".")[end]

fileName = ARGS[1]
filePath = dirname(ARGS[1])

# subtract extension and last dot from filename
fileNameBase = basename(ARGS[1])[1:end-(length(fileNameExtension)+1)]


nCells, nCellsPerBlock, nBlocks, nodeCoordinates, cubeIndices, numberDensity, minSize, maxSize, varNames = load_AMPS_data(fileName)
println(" - building nodes...")
nodes = build_nodes(nCells, nodeCoordinates, cubeIndices)
println(" - write result to .h5")

if isfile(joinpath(filePath, fileNameBase * ".h5"))
  println(" - file already exists, replacing it.")
  rm(joinpath(filePath, fileNameBase * ".h5"))
end

h5write(joinpath(filePath, fileNameBase * ".h5"), "oct/nCells", nCells)
h5write(joinpath(filePath, fileNameBase * ".h5"), "oct/nBlocks", nBlocks)
h5write(joinpath(filePath, fileNameBase * ".h5"), "oct/nCellsPerBlock", nCellsPerBlock)
h5write(joinpath(filePath, fileNameBase * ".h5"), "oct/nodes", nodes)
h5write(joinpath(filePath, fileNameBase * ".h5"), "oct/nodeCoordinates", nodeCoordinates)
h5write(joinpath(filePath, fileNameBase * ".h5"), "oct/cubeIndices", cubeIndices)
h5write(joinpath(filePath, fileNameBase * ".h5"), "oct/numberDensity", numberDensity)
h5write(joinpath(filePath, fileNameBase * ".h5"), "oct/minDustSize", minSize)
h5write(joinpath(filePath, fileNameBase * ".h5"), "oct/maxDustSize", maxSize)
h5write(joinpath(filePath, fileNameBase * ".h5"), "oct/varNames", varNames)
