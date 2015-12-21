using HDF5, JLD
using Types

initChildren = Array(Block, 8)
initCells = Array(Cell,1)

function build_octree(filePath, fileNameBase)
  nCells = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/nCells")
  nBlocks = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/nBlocks")
  nCellsPerBlock = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/nCellsPerBlock")
  nodes = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/nodes")
  nodeCoordinates = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/nodeCoordinates")
  cubeIndices = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/cubeIndices")
  numberDensity = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/numberDensity")
  varNames = h5read(joinpath(filePath, fileNameBase * ".h5"), "oct/varNames")

  xMax = maximum(nodes[:,:,1])
  yMax = maximum(nodes[:,:,2])
  zMax = maximum(nodes[:,:,3])

  halfSize = [xMax, yMax, zMax]
  root = [0.0, 0.0, 0.0]

  cellList = build_cells(nodes, cubeIndices, numberDensity, nCells)
  blocks = build_blocks(nBlocks, cellList, nodes, nCellsPerBlock)
  octree = Block(root, halfSize,1, initChildren, initCells, 5, 5, 5)
  if myid() == 1
    println(" - populating octree")
  end
  populate_octree(octree, blocks, nBlocks)

  nVars = size(numberDensity,1)

  return octree, nVars, varNames[4:end]
end


function is_out_of_bounds(oct, r)
  for i=1:3
    if ((r[i] > (oct.origin[i] + oct.halfSize[i])) | (r[i] < (oct.origin[i]-oct.halfSize[i])))
      return true
    end
  end
  return false
end

function octant_containing_point(point::Array{Float64,1}, block::Block)
  if !is_out_of_bounds(block, point)
    octant::Int64 = 1
    if (point[1] >= block.origin[1])
      octant += 4
    end
    if (point[2] >= block.origin[2])
      octant += 2
    end
    if (point[3] >= block.origin[3])
      octant += 1
    end
    return octant
  else
    return -1
  end
end

function block_containing_point(block::Block, point::Array{Float64,1})
  if !block.isLeaf
    oct = octant_containing_point(point, block)
    if oct == -1
      return false, block
    end
    block_containing_point(block.children[oct], point)
  else
    if !is_out_of_bounds(block, point)
      return true, block
    else
      return false, block
    end
  end
end

function cell_containing_point(oct::Block, point::Array{Float64, 1})
  foundBlock, block = block_containing_point(oct, point)
  if foundBlock
    nx = block.nx
    ny = block.ny
    nz = block.nz
    x = point[1] - block.cells[1].nodes[1,1]
    y = point[2] - block.cells[1].nodes[1,2]
    z = point[3] - block.cells[1].nodes[1,3]

    lx = block.halfSize[1] * 2.0 / nx
    ly = block.halfSize[2] * 2.0 / ny
    lz = block.halfSize[3] * 2.0 / nz

    fx = fld(x, lx)
    fy = fld(y, ly)
    fz = fld(z, lz)

    if fx > (nx-1.0)
        fx = nx - 1.0
    end
    if fy > (ny-1.0)
        fy = ny - 1.0
    end
    if fz > (nz-1.0)
        fz = nz - 1.0
    end

    cellIndex = round(Int, 1 + fx + fy*nx + fz*nx*ny)

    return true, block.cells[cellIndex]
  else
    return false, block.cells[1]
  end

end

function assign_triangles!(oct, allTriangles)
 point = zeros(Float64, 3)
 for tri in allTriangles
     for k=1:3
       for i=1:3
         point[i] = tri.nodes[i,k]
       end
       foundCell, cell = cell_containing_point(oct, point)
       if foundCell
         push!(cell.triangles, tri)
         cell.hasTriangles = true
         break
       end
     end
 end
end



function build_nodes(nCells::Int64, nodeCoordinates::Array{Float64,2},
                     cubeIndices::Array{Int64,2})
  nodes = Array(Float64, (nCells,8,3))

  for i=1:nCells
    for j=1:8
       for k=1:3
         @inbounds nodes[i,j,k] = nodeCoordinates[cubeIndices[i,j],k]
       end
    end
  end
  return nodes
end

function build_cells(nodes::Array{Float64,3}, cubeIndices::Array{Int64,2},
                     numberDensity::Array{Float64,2}, nCells::Int64)
  nVars = size(numberDensity, 1)
  nodeDensities = Array(Float64, (nVars, nCells, 8))
  origin = Array(Float64, (nCells,3))
  lCell = Array(Float64, 3)
  halfSize = Array(Float64, 3)
  allCells = Array(Cell, nCells)
  cellHasTriangles = false

  const lx = 0.0
  const ly = 0.0
  const lz = 0.0
  for i = 1:nCells
    for j = 1:8
      for k=1:nVars
        @inbounds nodeDensities[k, i, j] = numberDensity[k, cubeIndices[i,j]]
      end
    end
    volume = lx * ly * lz
    for k=1:3
      @inbounds lCell[k] = maximum(nodes[i,1:8,k]) - minimum(nodes[i,1:8,k])
      @inbounds halfSize[k] = lCell[k] / 2.0
      @inbounds origin[i,k] = nodes[i,1,k] + halfSize[k]
    end
    volume = lCell[1] * lCell[2] * lCell[3]
    @inbounds allCells[i] = Cell(vec(origin[i,1:3]),halfSize,
                                 reshape(nodes[i,1:8,1:3],(8,3)),
                                 volume,
                                 reshape(nodeDensities[1:nVars,i,1:8],(nVars, 8)),
                                 cellHasTriangles,
                                 Triangle[], 0)
  end
  return allCells
end

function build_blocks(nBlocks::Int64, allCells::Array{Cell,1}, nodes::Array{Float64,3},
                      nCellsPerBlock)
  lBlock = Array(Float64, (nBlocks,3))
  origin = Array(Float64, (nBlocks,3))
  blocks = Array(Block, nBlocks)

  for i=1:nBlocks
    blockCells = allCells[(i-1)*nCellsPerBlock+1:i*nCellsPerBlock]
    blockCellCoords = nodes[(i-1)*nCellsPerBlock+1: i*nCellsPerBlock, 1:8, 1:3]
    for j=1:3
      lBlock[i,j] = maximum(blockCellCoords[1:nCellsPerBlock,1:8,j]) - minimum(blockCellCoords[1:nCellsPerBlock,1:8,j])
      origin[i,j] = minimum(blockCellCoords[1:nCellsPerBlock,1:8,j]) + lBlock[i,j] / 2
    end
    @inbounds blocks[i] = Block(vec(origin[i,1:3]), vec(lBlock[i,1:3]).*0.5, 1, Array(Block, 8), blockCells,5,5,5)
  end
  return blocks
end

function splitBlock(node::Block)
  nx = node.nx
  ny = node.ny
  nz = node.nz

  xc1 = [node.origin[1] - node.halfSize[1]/2.0, node.origin[2] - node.halfSize[2]/2.0,  node.origin[3] - node.halfSize[3]/2.0]
  xc2 = [node.origin[1] - node.halfSize[1]/2.0, node.origin[2] - node.halfSize[2]/2.0,  node.origin[3] + node.halfSize[3]/2.0]
  xc3 = [node.origin[1] - node.halfSize[1]/2.0, node.origin[2] + node.halfSize[2]/2.0,  node.origin[3] - node.halfSize[3]/2.0]
  xc4 = [node.origin[1] - node.halfSize[1]/2.0, node.origin[2] + node.halfSize[2]/2.0,  node.origin[3] + node.halfSize[3]/2.0]
  xc5 = [node.origin[1] + node.halfSize[1]/2.0, node.origin[2] - node.halfSize[2]/2.0,  node.origin[3] - node.halfSize[3]/2.0]
  xc6 = [node.origin[1] + node.halfSize[1]/2.0, node.origin[2] - node.halfSize[2]/2.0,  node.origin[3] + node.halfSize[3]/2.0]
  xc7 = [node.origin[1] + node.halfSize[1]/2.0, node.origin[2] + node.halfSize[2]/2.0,  node.origin[3] - node.halfSize[3]/2.0]
  xc8 = [node.origin[1] + node.halfSize[1]/2.0, node.origin[2] + node.halfSize[2]/2.0,  node.origin[3] + node.halfSize[3]/2.0]

  node.children[1] = Block(xc1, node.halfSize/2, 1, Array(Block, 8), Array(Cell,1), nx, ny, nz)
  node.children[2] = Block(xc2, node.halfSize/2, 1, Array(Block, 8), Array(Cell,1), nx, ny, nz)
  node.children[3] = Block(xc3, node.halfSize/2, 1, Array(Block, 8), Array(Cell,1), nx, ny, nz)
  node.children[4] = Block(xc4, node.halfSize/2, 1, Array(Block, 8), Array(Cell,1), nx, ny, nz)
  node.children[5] = Block(xc5, node.halfSize/2, 1, Array(Block, 8), Array(Cell,1), nx, ny, nz)
  node.children[6] = Block(xc6, node.halfSize/2, 1, Array(Block, 8), Array(Cell,1), nx, ny, nz)
  node.children[7] = Block(xc7, node.halfSize/2, 1, Array(Block, 8), Array(Cell,1), nx, ny, nz)
  node.children[8] = Block(xc8, node.halfSize/2, 1, Array(Block, 8), Array(Cell,1), nx, ny, nz)
end

function getOctantContainingPoint(point::Array{Float64,1}, block::Block)
    octant::Int64 = 1
    if (point[1] >= block.origin[1])
      octant += 4
    end
    if (point[2] >= block.origin[2])
      octant += 2
    end
    if (point[3] >= block.origin[3])
      octant += 1
    end
  return octant
end

function insertChild(point::Array{Float64,1}, parent::Block, child::Block)
  if (parent.isLeaf && point != parent.origin)
    parent.isLeaf = false
    splitBlock(parent)
    octant = getOctantContainingPoint(point, parent)
    insertChild(point, parent.children[octant], child)
  elseif !parent.isLeaf
    octant = getOctantContainingPoint(point,  parent)
    insertChild(point, parent.children[octant], child)
  else
    parent.cells = child.cells
  end
end


function populate_octree(oct::Block, blocks::Array{Block, 1}, nBlocks::Int64)
  for i=1:nBlocks
    insertChild(blocks[i].origin, oct, blocks[i])
  end
  return oct
end


function findBlockContainingPoint(point::Array{Float64,1}, block::Block)
  if !block.isLeaf
    oct = getOctantContainingPoint(point, block)
    findBlockContainingPoint(point, block.children[oct])
  elseif block.isLeaf
    return block
  end
end

function findCellInBlock(block::Block, point::Array{Float64, 1})
  nx = 5.0
  ny = 5.0
  nz = 5.0
  x = point[1] - block.cells[1].nodes[1,1]
  y = point[2] - block.cells[1].nodes[1,2]
  z = point[3] - block.cells[1].nodes[1,3]
  lx = block.halfSize[1] * 2.0 / nx
  ly = block.halfSize[2] * 2.0 / ny
  lz = block.halfSize[3] * 2.0 / nz
  fx = fld(x, lx)
  fy = fld(y, ly)
  fz = fld(z, lz)

  if fx > (nx-1.0)
      fx = nx - 1.0
  end
  if fy > (ny-1.0)
      fy = ny-1.0
  end
  if fz > (nz-1.0)
      fz = nz-1.0
  end

  cellIndex = 1 + fx + fy*nx + fz*nx*ny
  return round(Int, cellIndex)
end

function get_factor(minSize, maxSize, r)
  f = Ref{Cdouble}(1.0)
  ccall((:ColumnIntegrationFactor, clib), Void, (Cdouble, Cdouble, Cdouble,
         Ptr{Cdouble}), minSize, maxSize, r, f)
  return f[]
end

function get_excitation_rate(coords)
  dist = sqrt(coords[1]*coords[1] + coords[2]*coords[2] + coords[3]*coords[3])
  exRate = 5.53e-11 * (100.64*1000 / dist)
  exRate = 5.53e-11
  return exRate
end

function triLinearInterpolation!(cell::Cell, point, data,
                                 minDustSize, maxDustSize, r)

  xd = (point[1] - cell.nodes[1,1]) / (cell.nodes[2,1] - cell.nodes[1,1])
  yd = (point[2] - cell.nodes[1,2]) / (cell.nodes[3,2] - cell.nodes[1,2])
  zd = (point[3] - cell.nodes[1,3]) / (cell.nodes[5,3] - cell.nodes[1,3])
  nVars = length(data)
  for i=1:nVars
    c00 = cell.densities[i,1] * (1-xd) + cell.densities[i, 2] * xd
    c10 = cell.densities[i,5] * (1-xd) + cell.densities[i, 6] * xd
    c01 = cell.densities[i,4] * (1-xd) + cell.densities[i, 3] * xd
    c11 = cell.densities[i,8] * (1-xd) + cell.densities[i, 7] * xd

    c0 = c00*(1-zd) + c10*zd
    c1 = c01*(1-zd) + c11*zd

    c = c0*(1-yd) + c1*yd

    if length(clib) > 1
      f = get_excitation_rate(point)
    else
      f = 1.0
    end
    data[i] = c * f
  end
end
