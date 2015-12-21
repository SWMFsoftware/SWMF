include("octree.jl")
include("io.jl")

function intersectBool(triangles, pStart, r)
  nTriangles = length(triangles)
  i = 0
  j = 0
  k = 0

  a = 0.0
  b = 0.0
  rI = 0.0

  dot_uv = 0.0
  dot_uu = 0.0
  dot_vv = 0.0
  dot_wv = 0.0
  dot_wu = 0.0

  divisor = 0.0
  sI = 0.0
  tI = 0.0

  # pI = intersect point
  pI = [0.,0.,0.]
  u = [0.,0.,0.]
  v = [0.,0.,0.]
  w = [0.,0.,0.]

  for i=1:nTriangles
    if dot(triangles[i].surfaceNormal, r) > 0.0
      continue
    end
    a = 0.0
    b = 0.0
    @simd for k=1:3
      @inbounds a = a + triangles[i].surfaceNormal[k] * (triangles[i].nodes[k,1] - pStart[k])
      @inbounds b = b + triangles[i].surfaceNormal[k] * r[k]
    end
    if ((a != 0.0) && (b != 0.0))
      rI = a / b
      if rI >= 0.0
        dot_uv = 0.0
        dot_uu = 0.0
        dot_vv = 0.0
        dot_wu = 0.0
        dot_wv = 0.0

        for k=1:3
          @inbounds pI[k] = pStart[k] + rI * r[k]
          @inbounds u[k] = triangles[i].nodes[k,2] - triangles[i].nodes[k,1]
          @inbounds v[k] = triangles[i].nodes[k,3] - triangles[i].nodes[k,1]
          @inbounds w[k] = pI[k] - triangles[i].nodes[k,1]

          @inbounds dot_uv = dot_uv + u[k]*v[k]
          @inbounds dot_uu = dot_uu + u[k]*u[k]
          @inbounds dot_vv = dot_vv + v[k]*v[k]
          @inbounds dot_wu = dot_wu + w[k]*u[k]
          @inbounds dot_wv = dot_wv + w[k]*v[k]
        end

        divisor = dot_uv * dot_uv - dot_uu * dot_vv
        sI = (dot_uv * dot_wv - dot_vv * dot_wu) / divisor
        tI = (dot_uv * dot_wu - dot_uu * dot_wv) / divisor

        if ((tI >= 0.0) && (sI >= 0.0) && (sI + tI < 1.0))
          return true
        end
      end
    end
  end
  return false
end

function iTriangleIntersect(triangles, pStart, r)
  nTriangles = length(triangles)
  iIntersectedTriangle = -1
  i = 0
  j = 0
  k = 0

  a = 0.0
  b = 0.0
  rI = 0.0

  dot_uv = 0.0
  dot_uu = 0.0
  dot_vv = 0.0
  dot_wv = 0.0
  dot_wu = 0.0

  divisor = 0.0
  sI = 0.0
  tI = 0.0

  # pI = intersect point
  pI = [0.,0.,0.]
  u = [0.,0.,0.]
  v = [0.,0.,0.]
  w = [0.,0.,0.]

  lIntersectMin = 0.0
  pIntersect = [0.0,0.0,0.0]

  for i=1:nTriangles
    if dot(triangles[i].surfaceNormal, r) > 0.0
      continue
    end
    a = 0.0
    b = 0.0
    @simd for k=1:3
      @inbounds a = a + triangles[i].surfaceNormal[k] * (triangles[i].nodes[k,1] - pStart[k])
      @inbounds b = b + triangles[i].surfaceNormal[k] * r[k]
    end
    if ((a != 0.0) && (b != 0.0))
      rI = a / b
      if rI >= 0.0
        dot_uv = 0.0
        dot_uu = 0.0
        dot_vv = 0.0
        dot_wu = 0.0
        dot_wv = 0.0

        for k=1:3
          @inbounds pI[k] = pStart[k] + rI * r[k]
          @inbounds u[k] = triangles[i].nodes[k,2] - triangles[i].nodes[k,1]
          @inbounds v[k] = triangles[i].nodes[k,3] - triangles[i].nodes[k,1]
          @inbounds w[k] = pI[k] - triangles[i].nodes[k,1]

          @inbounds dot_uv = dot_uv + u[k]*v[k]
          @inbounds dot_uu = dot_uu + u[k]*u[k]
          @inbounds dot_vv = dot_vv + v[k]*v[k]
          @inbounds dot_wu = dot_wu + w[k]*u[k]
          @inbounds dot_wv = dot_wv + w[k]*v[k]
        end
        lIntersect = 0.0
        for kk = 1:3
          lIntersect += (pI[kk] - pStart[kk]) * (pI[kk] - pStart[kk])
        end
        lIntersect = sqrt(lIntersect)

        divisor = dot_uv * dot_uv - dot_uu * dot_vv
        sI = (dot_uv * dot_wv - dot_vv * dot_wu) / divisor
        tI = (dot_uv * dot_wu - dot_uu * dot_wv) / divisor

        if ((tI >= 0.0) && (sI >= 0.0) && (sI + tI < 1.0))
          if (lIntersectMin == 0.0) | (lIntersect < lIntersectMin)
            lIntersectMin = lIntersect
            iIntersectedTriangle = i
            for kk =1:3
              pIntersect[kk] = pI[kk]
            end
          end
          continue
        end
      end
    end
  end
  return iIntersectedTriangle, pIntersect
end

function traverse_domain(oct::Block, r, r_hat, dr)
  foundCell, cell = cell_containing_point(oct, r)
  if foundCell
    dummyCell = deepcopy(cell)
    while true
      if is_out_of_bounds(oct, r)
        return -1, [0.0, 0.0, 0.0]
      end
      foundCell, cell = cell_containing_point(oct, r)
      if cell.hasTriangles
         iTriangle, pIntersect = iTriangleIntersect(cell.triangles, r, r_hat)
         if iTriangle != -1
           return iTriangle, pIntersect
         end
      end

      for k=1:3
        dr[k] = cell.halfSize[k] * r_hat[k]
        r[k] = r[k] + dr[k]
      end
    end
  else
    return -1, [0.0, 0.0, 0.0]
  end
end

function updateVectors!(r, r_hat, rStart, rPointing, i)
  for k=1:3
    @inbounds r[k] = rStart[k]
    @inbounds r_hat[k] = rPointing[k,i]
  end
end

function get_lMax(lMax, iTriangle, r, lIntersect, llMax)
  if iTriangle != -1
    for jj = 1:3
      @inbounds lMax += (r[jj] - lIntersect[jj]) * (r[jj] - lIntersect[jj])
    end
    lMax = sqrt(lMax)
  else
    lMax = llMax
  end
  return lMax
end

function reset_data!(nVars, data, dataNew, dataOld)
  for kk=1:nVars
    @inbounds data[kk] = 0.0
    @inbounds dataNew[kk] = 0.0
    @inbounds dataOld[kk] = 0.0
  end
end

function checkIfInShadow(doCheckShadow, cell, allTrianglesShadow,
                         r, r_hat_sun)
    if doCheckShadow == true
      if cell.isInShadow == 0
        isInShadow = intersectBool(allTrianglesShadow, r, r_hat_sun)
        if isInShadow
          cell.isInShadow = 1
        else
          cell.isInShadow = -1
        end
      elseif cell.isInShadow == 1
        isInShadow = true
      elseif cell.isInShadow == -1
        isInShadow = false
      end
    else
      isInShadow = false
    end
    return isInShadow
end


function update_r!(r, r_hat, dr, l)
  for k=1:3
    @inbounds dr[k] = r_hat[k] * l
    @inbounds r[k] = r[k] + dr[k]
  end
end

function norm_vec(r)
  sqrt(r[1]*r[1] + r[2]*r[2] + r[3]*r[3])
end

function doIntegration(oct::Block, rPointing, rStart, nVars, allTriangles,
                       doCheckShadow)
    if doCheckShadow
      meshFile = parseUserFile("meshFileShadow:")
      if length(meshFile) < 1
        println(" - shadow mesh file not defined, using regular mesh file")
        meshFile = parseUserFile("meshFile:")
      end
      nTrianglesShadow, allTrianglesShadow, totalSurfaceAreaShadow = load_ply_file(meshFile)
    end

    dataFileName = parseUserFile("dataFile:")
    minDustSize = h5read(dataFileName, "oct/minDustSize")
    maxDustSize = h5read(dataFileName, "oct/maxDustSize")

    nRays = size(rPointing, 2)
    columnDensity = 0.0
    distance = 0.0

    l = 0.0
    dr = Array(Float64, 3)
    r = Array(Float64, 3)
    r_hat = Array(Float64, 3)
    dr = zeros(Float64, 3)
    data = zeros(Float64, nVars)
    dataNew = zeros(Float64, nVars)
    dataOld = zeros(Float64, nVars)
    ccd = zeros(nVars, nRays)
    r_hat_sun = zeros(Float64, 3)
    rSun = zeros(Float64, 3)
    lMax = 0.0
    llMax = norm(oct.halfSize)*2.0
    distFromStart = 0.0
    mask = ones(Int64, nRays)
    userStr = lowercase(parseUserFile("pltBlankBody:"))
    if contains(userStr, "true") || contains(userStr, "yes")
      doBlankBody = true
    else
      doBlankBody = false
    end

    for i=1:nRays
      updateVectors!(r, r_hat, rStart, rPointing, i)
      iTriangle, lIntersect = iTriangleIntersect(allTriangles, r, r_hat)

      for k=1:3
        @inbounds r[k] = rStart[k]
      end
      #distance = sqrt(r[1]*r[1] + r[2]*r[2] + r[3]*r[3])
      distance = norm_vec(r)
      lMax = get_lMax(lMax, iTriangle, r, lIntersect, llMax)
      reset_data!(nVars, data, dataNew, dataOld)

      distFromStart = 0.0
      norm_rhat_sun = 0.0
      while (!is_out_of_bounds(oct, r) && (distFromStart <= lMax))
        if doCheckShadow
          spkpos!("SUN", et, "67P/C-G_CK", "NONE", "CHURYUMOV-GERASIMENKO", rSun)
        end
        for k=1:3
          @inbounds distFromStart += ((rStart[k] - r[k]) * (rStart[k] - r[k]))
          if doCheckShadow
            @inbounds r_hat_sun[k] = rSun[k] * 1000.0 - r[k]
            @inbounds norm_rhat_sun += r_hat_sun[k] * r_hat_sun[k]
          end
        end
        distFromStart = sqrt(distFromStart)

        if doCheckShadow
          norm_rhat_sun = sqrt(norm_rhat_sun)
          for k=1:3
            @inbounds r_hat_sun[k] = r_hat_sun[k] / norm_rhat_sun
          end
        end

        foundCell, cell = cell_containing_point(oct, r)
        l = norm_vec(r) / 30.0
        if doCheckShadow
          isInShadow = checkIfInShadow(doCheckShadow, cell, allTrianglesShadow,
                                      r, r_hat_sun)
        else
          isInShadow = false
        end

        if !isInShadow
          for k=1:nVars
            dataOld[k] = dataNew[k]
          end
          triLinearInterpolation!(cell, r, dataNew, minDustSize, maxDustSize,
                                  distFromStart)
          for k=1:nVars
            data[k] += ((dataOld[k] + dataNew[k]) / 2.0 * l)
          end
        end
        update_r!(r, r_hat, dr, l)
        distance = norm_vec(r)
      end

      if doBlankBody
        if iTriangle != -1
          mask[i] = 0
        end
      end
      for k=1:nVars
        ccd[k, i] = data[k]
      end
    end
    return ccd, mask
end

function doIntegrationParallel(oct::Block, rPointing, rStart, nVars, allTriangles,
                       doCheckShadow,iRays, et)
    if doCheckShadow
      meshFile = parseUserFile("meshFileShadow:")
      if length(meshFile) < 1
        println(" - shadow mesh file not defined, using regular mesh file")
        meshFile = parseUserFile("meshFile:")
      end
      nTrianglesShadow, allTrianglesShadow, totalSurfaceAreaShadow = load_ply_file(meshFile)
    end

    dataFileName = parseUserFile("dataFile:")
    minDustSize = h5read(dataFileName, "oct/minDustSize")
    maxDustSize = h5read(dataFileName, "oct/maxDustSize")

    nRays = size(rPointing, 2)
    columnDensity = 0.0
    distance = 0.0

    l = 0.0
    dr = Array(Float64, 3)
    r = Array(Float64, 3)
    r_hat = Array(Float64, 3)
    dr = zeros(Float64, 3)
    data = zeros(Float64, nVars)
    dataNew = zeros(Float64, nVars)
    dataOld = zeros(Float64, nVars)
    ccd = zeros(nVars, nRays)
    r_hat_sun = zeros(Float64, 3)
    rSun = zeros(Float64, 3)
    lMax = 0.0
    llMax = norm(oct.halfSize)*2.0
    distFromStart = 0.0
    mask = ones(Int64, nRays)
    userStr = lowercase(parseUserFile("pltBlankBody:"))
    if contains(userStr, "true") || contains(userStr, "yes")
      doBlankBody = true
    else
      doBlankBody = false
    end

    for i in iRays
      updateVectors!(r, r_hat, rStart, rPointing, i)
      iTriangle, lIntersect = iTriangleIntersect(allTriangles, r, r_hat)

      for k=1:3
        @inbounds r[k] = rStart[k]
      end
      #distance = sqrt(r[1]*r[1] + r[2]*r[2] + r[3]*r[3])
      distance = norm_vec(r)
      lMax = get_lMax(lMax, iTriangle, r, lIntersect, llMax)
      reset_data!(nVars, data, dataNew, dataOld)

      distFromStart = 0.0
      norm_rhat_sun = 0.0
      while (!is_out_of_bounds(oct, r) && (distFromStart <= lMax))
        if doCheckShadow
          spkpos!("SUN", et, "67P/C-G_CK", "NONE", "CHURYUMOV-GERASIMENKO", rSun)
        end
        for k=1:3
          @inbounds distFromStart += ((rStart[k] - r[k]) * (rStart[k] - r[k]))
          if doCheckShadow
            @inbounds r_hat_sun[k] = rSun[k] * 1000.0 - r[k]
            @inbounds norm_rhat_sun += r_hat_sun[k] * r_hat_sun[k]
          end
        end
        distFromStart = sqrt(distFromStart)

        if doCheckShadow
          norm_rhat_sun = sqrt(norm_rhat_sun)
          for k=1:3
            @inbounds r_hat_sun[k] = r_hat_sun[k] / norm_rhat_sun
          end
        end

        foundCell, cell = cell_containing_point(oct, r)
        l = norm_vec(r) / 30.0
        if doCheckShadow
          isInShadow = checkIfInShadow(doCheckShadow, cell, allTrianglesShadow,
                                      r, r_hat_sun)
        else
          isInShadow = false
        end

        if !isInShadow
          for k=1:nVars
            dataOld[k] = dataNew[k]
          end
          triLinearInterpolation!(cell, r, dataNew, minDustSize, maxDustSize,
                                  distFromStart)
          for k=1:nVars
            data[k] += ((dataOld[k] + dataNew[k]) / 2.0 * l)
          end
        end
        update_r!(r, r_hat, dr, l)
        distance = norm_vec(r)
      end

      if doBlankBody
        if iTriangle != -1
          mask[i] = 0
        end
      end
      for k=1:nVars
        ccd[k, i] = data[k]
      end
    end
    return ccd, mask
end

function dodo(iRays)
  doIntegrationParallel(oct, rPointing, rStart, nVars, allTriangles,
                         doCheckShadow_bool,iRays, et)
  end
