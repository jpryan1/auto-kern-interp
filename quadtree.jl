using RegionTrees
using StaticArrays: SVector
using Plots

import RegionTrees.Cell


include("structs.jl")

function node_subdivide(qt, node)

  split!(node, [ NodeData(), NodeData(), NodeData(), NodeData()])
  for i=1:2, j=1:2
    node[i,j].data.level = node.data.level+1
  end
  node_verts = vertices(node)
  mid = ( node_verts[2,2] - node_verts[1,1] ) / 2. + node_verts[1,1]
  
  if node.data.level + 1 > length(qt.levels)
    new_level = QuadTreeLevel([node[1,1], node[1,2], node[2,1], node[2,2]])
    push!(qt.levels, new_level)
  else
    for i=1:2, j=1:2
      push!(qt.levels[node.data.level+1].nodes , node[i,j])
    end
  end
  
  
  for index in 1:qt.solution_dimension:length(node.data.dof_lists.original_box)
  
    matrix_index = node.data.dof_lists.original_box[index]
  
    points_vec_index = 2*matrix_index - 1
    if qt.solution_dimension == 2
      points_vec_index = matrix_index
    end
    
    x = qt.points[points_vec_index]
    y = qt.points[points_vec_index + 1]

    if x <= mid[1] && y <= mid[2]

      for i in 1:qt.solution_dimension
        push!(node[1,1].data.dof_lists.original_box, matrix_index + i - 1)
      end
    elseif x <= mid[1] && y > mid[2]

      for i in 1:qt.solution_dimension
        push!(node[1,2].data.dof_lists.original_box, matrix_index + i - 1)
      end
    elseif x > mid[1] && y <= mid[2]

      for i in 1:qt.solution_dimension
        push!(node[2,1].data.dof_lists.original_box, matrix_index + i - 1)
      end
    else

      for i in 1:qt.solution_dimension
        push!(node[2,2].data.dof_lists.original_box, matrix_index + i - 1)
      end
    end
  end
  
  for i in 1:2, j in 1:2
    if length(node[i,j].data.dof_lists.original_box) > 32
      node_subdivide(qt, node[i, j])
    end
  end
end


function add_point(qt, node, point, point_num)
  if qt.solution_dimension == 1
    push!(node.data.dof_lists.original_box, point_num)
  elseif qt.solution_dimension == 2
    push!(node.data.dof_lists.original_box, 2*(point_num-1) + 1)
    push!(node.data.dof_lists.original_box, 2*(point_num-1) + 2)
  end
  
  node_verts = vertices(node)
  mid = ( node_verts[2,2] - node_verts[1,1] ) / 2. + node_verts[1,1]

  if isleaf(node)
    if length(node.data.dof_lists.original_box) > 32  # MAX LEAF DOFS
      node_subdivide(qt, node)
    end
  else
    if point[1] <= mid[1] && point[2] <= mid[2]
      add_point(qt, node[1,1], point, point_num)
    elseif point[1] > mid[1] && point[2] <= mid[2]
      add_point(qt, node[2,1], point, point_num)
    elseif point[1] <= mid[1] && point[2] > mid[2]
      add_point(qt, node[1,2], point, point_num)
    else
      add_point(qt, node[2,2], point, point_num)
    end
  end
end


function compute_neighbors(qt)
  for level in 2:length(qt.levels)
    current_level = qt.levels[level]
    for node_a in current_level.nodes
      for sibling in node_a.parent.children
        if sibling != node_a
          push!(node_a.data.neighbors, sibling)
        end
      end
      for parents_neighbor in node_a.parent.data.neighbors
        if isleaf(parents_neighbor)
          continue
        end
        for cousin in parents_neighbor.children
          if(cousin.data.level != node_a.data.level)
            continue
          end
          node_a_verts = vertices(node_a)
          cousin_verts = vertices(cousin)
          dist = norm(node_a_verts[1,1] - cousin_verts[1,1])
          side_length = norm(node_a_verts[1,1] - node_a_verts[1,2])
          if(dist < side_length * sqrt(2) + 1e-6)
            push!(node_a.data.neighbors, cousin)
          end
        end
      end
      if isleaf(node_a)
        for neighbor in node_a.data.neighbors
          if neighbor.data.level != node_a.data.level
            continue
          end
          if isleaf(neighbor)
            continue
          end
          for child in neighbor.children
            get_descendent_neighbors(qt, node_a, child)
          end
        end
      end
    end
  end
end


function get_descendent_neighbors(qt, big, small)

  big_verts = vertices(big)
  big_top = big_verts[2,2][2]
  big_bottom = big_verts[1,1][2]
  big_left = big_verts[1,1][1]
  big_right = big_verts[2,2][1]

  for small_vert in vertices(small)
    if small_vert[1] > big_left && small_vert[1] < big_right
      if abs(small_vert[2] - big_top) < 1e-14
        push!(big.data.neighbors, small)
        push!(small.data.neighbors, big)
        break
      elseif abs(small_vert[2] - big_bottom) < 1e-14
        push!(big.data.neighbors, small)
        push!(small.data.neighbors, big)
        break
      end
    elseif small_vert[2] < big_top && small_vert[2] > big_bottom
      if abs(small_vert[1] - big_left) < 1e-14
        push!(big.data.neighbors, small)
        push!(small.data.neighbors, big)
        break
      elseif abs(small_vert[1] - big_right) < 1e-14
        push!(big.data.neighbors, small)
        push!(small.data.neighbors, big)
        break
      end
    end
  end
  if !isleaf(small)
    for child in small.children
      get_descendent_neighbors(qt, big, child)
    end
  end
end


function initialize_tree(points)
  point_min = -0.01 + minimum(points)
  point_max = 0.01 + maximum(points)
  root_data = NodeData()
  root = Cell(SVector(point_min,point_min), SVector(point_max-point_min, point_max-point_min), root_data)
  first_level = QuadTreeLevel([root])
  levels = [first_level]


  qt = QuadTree(root, levels, points, Matrix{Float64}(undef, 0, 0), 1)
  for point_num in 1:(length(points)÷2)
    add_point(qt, root, [points[2*point_num-1], points[2*point_num]], point_num)
  end
  compute_neighbors(qt)

  plt = plot(xlim=(point_min, point_max), ylim=(point_min, point_max), legend=nothing)
  for leaf in allleaves(root)
      v = hcat(collect(vertices(leaf.boundary))...)
      plot!(plt, v[1,[1,2,4,3,1]], v[2,[1,2,4,3,1]])
  end
  gui()
  return qt
end
