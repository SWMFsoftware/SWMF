using Types
module Triangles

export build_triangles,
       calculate_surface_normals,
       calculate_tri_centers,
       calculate_tri_areas,
       assign_triangles!



function calculate_surface_normals(nodeCoords, triIndices, nTriangles)

  n_hat = zeros(Float64, 3, nTriangles)
  vi = zeros(Float64, 3)
  vj = zeros(Float64, 3)
  vk = zeros(Float64, 3)

  for ii=1:nTriangles
    i = triIndices[1, ii]
    j = triIndices[2, ii]
    k = triIndices[3, ii]

    vi = vec(nodeCoords[1:3, i])
    vj = vec(nodeCoords[1:3, j])
    vk = vec(nodeCoords[1:3, k])
    r = cross(vj-vi, vk-vi)
    r = r/norm(r)
    for kk = 1:3
      n_hat[kk, ii] = r[kk]
    end
  end

  return n_hat

end

function build_triangles(nodeCoords, triIndices, nTriangles)
  triangles = zeros(Float64, 3, 3, nTriangles)
  for i=1:nTriangles
    for j=1:3
      for k=1:3
        triangles[k,j,i] = nodeCoords[k,triIndices[j,i]]
      end
    end
  end
  return triangles
end

function calculate_tri_centers(triangles, nTriangles)
  triCenters = zeros(Float64, 3, nTriangles)
  for i=1:nTriangles
    for j=1:3
      triCenters[j,i] = sum(triangles[j,1:3,i])/3.0
    end
  end
  return triCenters
end

function calculate_tri_areas(triangles, nTriangles)
  triAreas = zeros(Float64, nTriangles)
  for i=1:nTriangles
    P = vec(triangles[1:3,2,i] - triangles[1:3,1,i])
    Q = vec(triangles[1:3,3,i] - triangles[1:3,1,i])
    S = sqrt(sum(cross(P,Q).^2))
    triAreas[i] = 0.5 * S
  end
  return triAreas
end

end
