nu1 = 50
nu2 = 100
u1grid = zeros(nu1)
u2grid = zeros(nu2)
u1i = 0.0
u1f = 0.5
u2i = -0.25
u2f = 0.75
for i=1:nu1
  u1grid[i] = (i-0.5)*(u1f-u1i)/nu1+u1i
end
for i=1:nu2
  u2grid[i] = (i-0.5)*(u2f-u2i)/nu2+u2i
end

uch1_grid = zeros(nu1,nu2)
uch2_grid = zeros(nu1,nu2)
for j=1:size(u11)[1]
  if u11[j] > u1i && u11[j] < u1f && u21[j] > u2i && u21[j] < u2f
    iu1 = ceil(Int64,(u11[j]-u1i)/(u1f-u1i)*nu1)
    iu2 = ceil(Int64,(u21[j]-u2i)/(u2f-u2i)*nu2)
    uch1_grid[iu1,iu2] +=1.0
  end
  if u12[j] > u1i && u12[j] < u1f && u22[j] > u2i && u22[j] < u2f
    iu1 = ceil(Int64,(u12[j]-u1i)/(u1f-u1i)*nu1)
    iu2 = ceil(Int64,(u22[j]-u2i)/(u2f-u2i)*nu2)
    uch2_grid[iu1,iu2] +=1.0
  end
end

# Set levels for each limb-darkening parameter:
conf = [0.95, 0.683]; nconf=2
uch1_grid ./= float(sum(uch1_grid))
uch2_grid ./= float(sum(uch2_grid))
uch1_grid_sort = sort(vec(uch1_grid))
uch2_grid_sort = sort(vec(uch2_grid))
i0 = nu1*nu2
uch1_grid_cum = uch1_grid_sort[i0]
level_uch1 = zeros(nconf)
while uch1_grid_cum < conf[1] && i0 > 1
  global i0 -= 1
  for k=1:nconf
    if uch1_grid_cum < conf[k] && (uch1_grid_cum+uch1_grid_sort[i0]) > conf[k]
      level_uch1[k] = 0.5*(uch1_grid_sort[i0]+uch1_grid_sort[i0+1])
    end
  end
  global uch1_grid_cum += uch1_grid_sort[i0]
end
uch2_grid_cum = uch2_grid_sort[i0]
level_uch2 = zeros(nconf)
i0 = nu1*nu2
while uch2_grid_cum < conf[1] && i0 > 1
  global i0 -= 1
  for k=1:nconf
    if uch2_grid_cum < conf[k] && (uch2_grid_cum+uch2_grid_sort[i0]) > conf[k]
      level_uch2[k] = 0.5*(uch2_grid_sort[i0]+uch2_grid_sort[i0+1])
    end
  end
  global uch2_grid_cum += uch2_grid_sort[i0]
end
cs = contour(u1grid,u2grid,uch1_grid',levels=level_uch1,colors="r",linewidth=2.0)
cs = contour(u1grid,u2grid,uch2_grid',levels=level_uch2,colors="g",linewidth=2.0)
