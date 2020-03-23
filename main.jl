module auto_kern_interp

using LinearAlgebra

include("quadtree.jl");
include("skel-interp.jl");

function kernel(x,y)
  if(x==y)
    return 2
  end
  return 1.0/norm(x-y)^2
end


function main() 

  num_points = 2^(8)
  points_vec = zeros(num_points*2)
  for idx in 1:num_points
    ang = 2.0*idx*pi/num_points
    
    points_vec[2*idx-1] = 0.375*cos(ang)*(sin(ang)+4)
    points_vec[2*idx] = 0.375*sin(ang)*(sin(ang)+4)
  end

  quadtree = initialize_tree(points_vec)
  skel_interp(kernel, quadtree)

  f = randn((num_points, 1))
  mu = randn((num_points, 1))
  solve(quadtree, mu, f)
 
  points_arr = [[points_vec[2*i-1], points_vec[2*i]] for i in 1:num_points]
  forward_mat = kernel.(points_arr, permutedims(points_arr)) 

  forward = forward_mat * mu
  println("Err ", norm(forward - f)/norm(f))
end

end

auto_kern_interp.main()
