using HDF5, JLD
using PyPlot
using Spice
using Instrument

include("io.jl")
include("octree.jl")
@everywhere include("raytrace.jl")

global const clib = parseUserFile("clibFile:")
@show(clib)

etStr = ARGS[1]
metaFile = parseUserFile("kernelFile:")
try
  furnsh(metaFile)
catch e
  println("spice error")
  println(e)
  exit()
end
et = str2et(etStr)
if length(ARGS) == 2
  instrumentName = ARGS[2]
  instrument = rosetta_instruments[instrumentName]
  nPixelsX = instrument.nPixelsX
  nPixelsY = instrument.nPixelsY
  phiX = instrument.phiX
  phiY = instrument.phiY
  rotMat = pxform(instrument.frame, "67P/C-G_CK", et)
  rotMat = rotMat'
else
  nPixelsX = parse(Int, ARGS[2])
  nPixelsY = parse(Int, ARGS[3])
  phiX = parse(Float64, ARGS[4])
  phiY = parse(Float64, ARGS[5])

  rotMat = eye(3)
end

const fileName = parseUserFile("dataFile:")
const meshFile = parseUserFile("meshFile:")
const doCheckShadow = parseUserFile("doCheckShadow:")
if contains(doCheckShadow, "y")
  doCheckShadow_bool = true
else
  doCheckShadow_bool = false
end

doBlankBody = false
userStr = lowercase(parseUserFile("pltBlankBody:"))
if contains(userStr, "true") || contains(userStr, "yes")
  doBlankBody = true
else
  doBlankBody = false
end

const filePath = dirname(fileName)
fileNameExtension = split(basename(fileName), ".")[end]
fileNameBase = basename(fileName)[1:end-(length(fileNameExtension)+1)]

################################################################################
# load data from spice
rRos_km, lt = spkpos("ROSETTA", et, "67P/C-G_CK", "NONE", "CHURYUMOV-GERASIMENKO")
println(" - observer distance from coordinate center : ", norm(rRos_km), " km")
rStart = rRos_km .* 1000

################################################################################
# create pointing vectors
lx = 2 * sin(phiX / 180.0 * pi)
ly = 2 * sin(phiY / 180.0 * pi)

rPointing = ones(Float64, 3, nPixelsX * nPixelsY)
z = 1.0
k = 1
for (j,y) = enumerate(linspace(-ly / 2.0, ly / 2.0, nPixelsY))
  for (i,x) = enumerate(linspace(-lx / 2.0, lx / 2.0, nPixelsX))
    norm_rPointing = sqrt(x*x + y*y + z*z)
    rPointing[1,k] = x / norm_rPointing
    rPointing[2,k] = y / norm_rPointing
    rPointing[3,k] = z / norm_rPointing
    k += 1
  end
end

# rotate pointing vectors to comet centric frame
for k=1:nPixelsX*nPixelsY
  rPointing[:,k] = rotMat * vec(rPointing[:,k])
end


oct, nVars = build_octree(filePath, fileNameBase)
println(" - nVars : ", nVars)
nTriangles, allTriangles, totalSurfaceArea = load_ply_file(meshFile)

assign_triangles!(oct, allTriangles)

@time ccd, mask = doIntegration(oct, rPointing, rStart, nVars, allTriangles,
                          doCheckShadow_bool)

ccd = reshape(ccd, nVars, nPixelsX, nPixelsY)
ccd_sum = zeros(nPixelsX, nPixelsY)
if doBlankBody
  mask = reshape(mask, nPixelsX, nPixelsY)
  border_mask = get_border(mask)
  custom_mask!(mask)
end

################################################################################
# check if user set plotting options
################################################################################
cmap = ColorMap("hot")
if length(parseUserFile("pltColorMap:")) > 0
  cmapUser = parseUserFile("pltColorMap:")
  try
    cmap = ColorMap(cmapUser)
  catch
    print_with_color(:red, " - your colormap was not found. using default 'hot'\n")
    println(" - some valid choices are: afmhot, autumn, bone, cool")
    println(" - copper, gist_heat, gray, pink, spring, summer, winter")
    println(" - Blues, Greens, Oranges, Reds, YlGn, BuPu, Greys, YlGnBu")
  end
end

nLevels = 32
if length(parseUserFile("pltLevels:")) > 0
  nLevels = parse(Int, parseUserFile("pltLevels:"))
end

pltTitle = "I love Julia (by Valeriy T.)"
if length(parseUserFile("pltTitle:")) > 0
  pltTitle = parseUserFile("pltTitle:")
end

fontSize = 12
if length(parseUserFile("pltFontSize:")) > 0
  fontSize = parse(Int, parseUserFile("pltFontSize:"))
end

for i=1:nVars
  figure()
  ccdPlt = reshape(ccd[i,:,:], nPixelsX, nPixelsY)
  contourf(log10(ccdPlt), nLevels, cmap=cmap)
  colorbar()
  if doBlankBody
    contourf(mask, levels=[-0.1, 0.1], colors=("w"))
    contourf(border_mask, levels=[0.9, 1.1], colors=("k"))
  end
  xlabel("Pixel number", size=fontSize)
  ylabel("Pixel number", size=fontSize)
  title(pltTitle, size=fontSize)
  for ix = 1:nPixelsX
    for iy = 1:nPixelsY
       ccd_sum[ix, iy] += ccdPlt[ix, iy]
     end
   end
end
figure()
contourf(log10(ccd_sum), nLevels, cmap=cmap)
xlabel("Pixel number", size=fontSize)
ylabel("Pixel number", size=fontSize)
title("Sum", size=fontSize)
colorbar()


show()

writedlm(joinpath(filePath, "ccd.dat"), ccd)
