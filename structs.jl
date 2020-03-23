
import RegionTrees.Cell

mutable struct DofLists
  original_box::Array{Int64}
  active_box::Array{Int64}
  skel::Array{Int64}
  skelnear::Array{Int64}
  near::Array{Int64}
  redundant::Array{Int64}
  perm::Array{Int64}
end
DofLists() = DofLists([],[],[],[],[], [], [])


mutable struct NodeData
  T::Matrix{Float64}
  L::Matrix{Float64}
  U::Matrix{Float64}
  X_rr::Matrix{Float64}
  schur_update::Matrix{Float64}
  level::Int64
  dof_lists::DofLists
  compressed::Bool
  neighbors::Array{Cell{NodeData, 2, Float64, 4}}
end
NodeData() = NodeData(Matrix{Float64}(undef, 0, 0),Matrix{Float64}(undef, 0, 0),Matrix{Float64}(undef, 0, 0),Matrix{Float64}(undef, 0, 0),Matrix{Float64}(undef, 0, 0), 1, DofLists(), false, [])


mutable struct QuadTreeLevel
  nodes::Array{Cell{NodeData, 2, Float64, 4}}
end


mutable struct QuadTree
  root::Cell{NodeData, 2, Float64, 4}
  levels::Array{QuadTreeLevel}
  points::Array{Float64}
  allskel_mat::Matrix{Float64}
  solution_dimension::Int64
end
