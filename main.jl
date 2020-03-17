

using LinearAlgebra

include("quadtree.jl");
include("skel-interp.jl");

function kernel(x,y)
  if(norm(x-y)<1e-6)
    return 2
  end
  return exp(norm(x-y)^2)
end


function main()

  x_points = [[x,y] for x in 0.1:0.01:0.9, y in 0.1:0.01:0.9]
  points_vec = zeros(16000)
  for ang in 1:8000
    points_vec[2*ang-1] = cos((ang)*pi/4000.)
    points_vec[2*ang] = sin((ang)*pi/4000.)
  end

  quadtree = initialize_tree(points_vec)
  skel_interp(kernel, quadtree)

  f = ones((8000, 1))
  mu = ones((8000, 1))
  solve(quadtree, mu, f)
  # mu = quadtree.allskel_mat \ f
  points_arr = [[points_vec[2*i-1], points_vec[2*i]] for i in 1:8000]
  forward_mat = kernel.(points_arr, permutedims(points_arr)) 

  forward = forward_mat * mu
  println("Err ", norm(forward - f))
  println("Done.")
end

main()
