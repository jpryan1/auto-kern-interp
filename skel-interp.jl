
function update_active_dofs(tree)
  for level in length(tree.levels):-1:1
    update_active_dofs(tree, level)
  end
end


function update_active_dofs(tree, level::Int32)
  current_level = tree.levels[level]
  for node in current_level.nodes
    if node.data.compressed
      continue
    end
    update_active_dofs(tree, node)
  end
  # Next, get all active near dofs from neighbors
  for node_a in current_level.nodes
    node_a.data.dof_lists.near = []
    for neighbor in node_a.data.neighbors
      if neighbor.data.level > node_a.data.level
        continue
      end
      if isleaf(neighbor)
        append!(node_a.data.dof_lists.near, neighbor.data.dof_lists.original_box)
      else
        append!(node_a.data.dof_lists.near, neighbor.data.dof_lists.active_box)
      end
    end
  end
end


function update_active_dofs(tree, node)
  node.data.dof_lists.skel = []
  node.data.dof_lists.skelnear = []
  node.data.dof_lists.redundant = []
  node.data.dof_lists.active_box = []

  if !isleaf(node)
    for child in node.children
      if child.data.compressed
        append!(node.data.dof_lists.active_box, child.data.dof_lists.skel)
      else
        append!(node.data.dof_lists.active_box, child.data.dof_lists.active_box)
      end
    end
  else
    node.data.dof_lists.active_box = copy(node.data.dof_lists.original_box)
  end
end


function get_all_schur_updates(indices, node)
  updates = zeros(length(indices), length(indices))
  if !isleaf(node) 
    get_descendents_updates(updates, indices, node)
  end
  return updates
end


function get_descendents_updates(updates, BN, node)
  for child in node.children
    if child.data.compressed
      get_update(updates, BN, child)
    end
    if !isleaf(child)
      get_descendents_updates(updates, BN, child)
    end
  end
end


function get_update(updates, BN, node)
  sn_indices = []
  BN_indices = []
  for sn_idx in 1:length(node.data.dof_lists.skelnear)
    for BN_idx in 1:length(BN)
      if BN[BN_idx] == node.data.dof_lists.skelnear[sn_idx]
        push!(sn_indices, sn_idx)
        push!(BN_indices, BN_idx)
      end
    end
  end
  num_shared_by_both = length(BN_indices)
  for i in 1:num_shared_by_both
    for j in 1:num_shared_by_both
      updates[BN_indices[i], BN_indices[j]] += node.data.schur_update[sn_indices[i], sn_indices[j]]
    end
  end
end

  
function id(A::Matrix{Float64})
# in place? Info return?
  Q, R, P = qr(A, Val(true))
  thresh = 1e-6* abs(R[1,1])
  best = 0
  for i in 2:length(P)
    if(abs(R[i,i]) < thresh)
      best = i
      break
    end
  end
  if best >= length(P)
    return nothing, nothing
  end
  Z = R[1:best, 1:best] \ R[1:best, (best+1):end]
  return Z, P
end


function get_id_mat(kernel, tree, node)
    #TODO Make this into a function
  active_box_points = [[tree.points[2*idx-1],tree.points[2*idx]] 
                        for idx in node.data.dof_lists.active_box]
  # If level one or two, use full offdiag, no proxy
  outside_box_points = []
  if node.data.level == 1
    outside_box_points = [[tree.points[2*idx-1],tree.points[2*idx]] 
                          for idx in node.data.dof_lists.near]
    id_mat = Matrix{Float64}(undef, 2*length(outside_box_points), length(active_box_points))
    
    id_mat[1:length(outside_box_points), :] = 
      kernel.(outside_box_points, permutedims(active_box_points))
    
    id_mat[(length(outside_box_points)+1):end, :] = 
      transpose(kernel.(active_box_points, permutedims(outside_box_points)))
    
    return id_mat
  end
  if node.data.level == 2
    for level_node in tree.levels[2].nodes
      if level_node == node
        continue
      end
      append!(outside_box_points, [[tree.points[2*idx-1],tree.points[2*idx]] 
              for idx in level_node.data.dof_lists.active_box])
    end
    for level_node in tree.levels[1].nodes
      if !isleaf(level_node)
        continue
      end
      append!(outside_box_points, [[tree.points[2*idx-1],tree.points[2*idx]] 
              for idx in level_node.data.dof_lists.active_box])
    end
    id_mat = Matrix{Float64}(undef, 2*length(outside_box_points), length(active_box_points))
    
    id_mat[1:length(outside_box_points), :] = 
      kernel.(outside_box_points, permutedims(active_box_points))
   
    id_mat[(length(outside_box_points)+1):end, :] = 
      transpose(kernel.(active_box_points, permutedims(outside_box_points)))
    
    return id_mat
  end


  node_verts = vertices(node)
  side_length = norm(node_verts[1,1] - node_verts[1,2])
  mid = ( node_verts[2,2] - node_verts[1,1] ) / 2. + node_verts[1,1]
  nearby_points = [[tree.points[2*idx-1],tree.points[2*idx]] for idx in node.data.dof_lists.near]
  for point in nearby_points
    if(norm(point-mid) < 1.5*side_length)
      push!(outside_box_points, point)
    end
  end
  pxy_points = [mid + 1.5*side_length*[cos(ang), sin(ang)]
                for ang in (2*pi/128.):(2*pi/128.):(2*pi)]

  idx_a = length(outside_box_points)
  idx_b = 2*length(outside_box_points)
  idx_c = 2*length(outside_box_points) + length(pxy_points)
  idx_d = 2*length(outside_box_points) + 2*length(pxy_points)

  # id_mat = Matrix{Float64}(0., idx_d, length(active_box_points))
  id_mat = zeros( idx_d, length(active_box_points))
  
  id_mat[1:idx_a, :] = 
    kernel.(outside_box_points, permutedims(active_box_points))
  # id_mat[(idx_a+1):idx_b, :] = 
  #   transpose(kernel.(active_box_points, permutedims(outside_box_points))) 

  # id_mat[(idx_b+1):idx_c, :] = 
  #   kernel.(pxy_points, permutedims(active_box_points))
  # id_mat[(idx_c+1):idx_d, :] = 
  #   transpose(kernel.(active_box_points, permutedims(pxy_points))) 
  return id_mat

end


function id_compress(kernel, tree, node)
  pxy = get_id_mat(kernel, tree, node)
  T, P = id(pxy);
  if P == nothing
    return 0
  end
  node.data.T = T
  numskel = size(node.data.T, 1)

  outside_box_points = []
  for level_node in tree.levels[node.data.level].nodes
    if level_node == node
      continue
    end
    append!(outside_box_points, [[tree.points[2*idx-1],tree.points[2*idx]] 
            for idx in level_node.data.dof_lists.active_box])
  end
  # active_box_points = [[tree.points[2*idx-1],tree.points[2*idx]] 
  #                       for idx in node.data.dof_lists.active_box]
 

  node.data.dof_lists.perm = P
  node.data.dof_lists.skel = node.data.dof_lists.active_box[P[1:numskel]]
  node.data.dof_lists.redundant = node.data.dof_lists.active_box[P[(numskel+1):end]]
  node.data.dof_lists.skelnear = copy(node.data.dof_lists.skel)
       

  redund_points =  [[tree.points[2*idx-1],tree.points[2*idx]] 
                        for idx in node.data.dof_lists.redundant]
  skel_points =  [[tree.points[2*idx-1],tree.points[2*idx]] 
                        for idx in node.data.dof_lists.skel]
 
  KfsT = kernel.(outside_box_points, permutedims(skel_points)) * T
  err = norm(KfsT - kernel.(outside_box_points, permutedims(redund_points))) / norm(kernel.(outside_box_points, permutedims(redund_points)))
  println("ID err ", err) 
     
  return numskel
end


function decouple(kernel, tree, node)
  num_skel = size(node.data.T, 1)
  num_redundant = size(node.data.T, 2)
  BN = copy(node.data.dof_lists.active_box)
  BN_points = [[tree.points[2*idx-1],tree.points[2*idx]] for idx in BN]
  schur_updates = get_all_schur_updates(BN, node)
  K_BN = kernel.(BN_points, permutedims(BN_points)) - schur_updates
  s = node.data.dof_lists.perm[1:num_skel]
  sn = node.data.dof_lists.perm[1:num_skel]
  r = node.data.dof_lists.perm[(num_skel+1):(num_skel+num_redundant)]
  K_BN_r_sn = K_BN[r, s] - transpose(node.data.T)*K_BN[s,s]
  node.data.X_rr = K_BN[r,r] - transpose(node.data.T)*K_BN[s,r] - K_BN_r_sn*node.data.T
  K_BN_sn_r = K_BN[s, r] - K_BN[s, s] * node.data.T
  node.data.L = transpose(transpose(node.data.X_rr) \ transpose(K_BN_sn_r)) 
  node.data.U = node.data.X_rr \ K_BN_r_sn
  node.data.schur_update = node.data.L * K_BN_r_sn
  node.data.compressed = true
end


function skel_interp(kernel, tree)
  lvls = length(tree.levels)
  for level in lvls:-1:2
    println("level ", level)
    update_active_dofs(tree, level)
    current_level = tree.levels[level]
    for current_node in current_level.nodes
      if current_node.data.compressed || length(current_node.data.dof_lists.active_box) < 16 #min dofs to compress
        continue
      end
      if id_compress(kernel, tree, current_node) == 0
        continue
      end
      decouple(kernel, tree, current_node)
    end
  end
  update_active_dofs(tree)
  allskel = tree.root.data.dof_lists.active_box
  
  if length(allskel) > 0
    allskel_updates = get_all_schur_updates(allskel, tree.root)
    allskel_points = [[tree.points[2*allskel[i]-1], tree.points[2*allskel[i]]] for i in 1:length(allskel)]
    tree.allskel_mat = kernel.(allskel_points, permutedims(allskel_points)) - allskel_updates
  end

  #   tree->allskel_mat.LU_factorize(&tree->allskel_mat_lu,
  #                                 &tree->allskel_mat_piv);

end


function solve(quadtree, x, b)
  x[:] = copy(b)
  lvls = length(quadtree.levels)
  all_nodes = []
  for level in lvls:-1:1
    current_level = quadtree.levels[level]
    append!(all_nodes, current_level.nodes)
  end
  for level in lvls:-1:1
    current_level = quadtree.levels[level]
    for current_node in current_level.nodes
      if !current_node.data.compressed
        continue
      end

      x[current_node.data.dof_lists.redundant] += (-transpose(current_node.data.T)) * x[current_node.data.dof_lists.skel]
      x[current_node.data.dof_lists.skelnear] += (-current_node.data.L) * x[current_node.data.dof_lists.redundant]
    end
  end
  for n in 1:length(all_nodes)
    current_node = all_nodes[n]
    if length(current_node.data.dof_lists.redundant) == 0 || !current_node.data.compressed
      continue
    end
    x[current_node.data.dof_lists.redundant] =  current_node.data.X_rr \ x[current_node.data.dof_lists.redundant]
  end
  allskel = quadtree.root.data.dof_lists.active_box
  if length(allskel) > 0
    x[allskel] = quadtree.allskel_mat \ x[allskel]
  end
  for level in 1:lvls
    current_level = quadtree.levels[level]
    for current_node in current_level.nodes
      if !current_node.data.compressed
        continue
      end
      x[current_node.data.dof_lists.redundant] += (-current_node.data.U) * x[current_node.data.dof_lists.skelnear]
      x[current_node.data.dof_lists.skel] += (-current_node.data.T) * x[current_node.data.dof_lists.redundant]
    end
  end
end
