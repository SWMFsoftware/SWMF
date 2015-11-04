module Types

export Triangle,
       Cell,
       Block

type Triangle
  id::Int64
  center::Array{Float64,1}
  nodes::Array{Float64,2}
  area::Float64
  surfaceNormal::Array{Float64,1}
end

type Cell
  origin::Array{Float64,1}
  halfSize::Array{Float64,1}
  nodes::Array{Float64,2}
  volume::Float64
  densities::Array{Float64, 2}
  hasTriangles::Bool
  triangles::Vector{Triangle}
  isInShadow::Int64
end

type Block
  origin::Array{Float64, 1}
  halfSize::Array{Float64, 1}
  isLeaf::Bool
  children::Array{Block, 1}
  cells::Array{Cell,1}
  nx::Int64
  ny::Int64
  nz::Int64
end

end
