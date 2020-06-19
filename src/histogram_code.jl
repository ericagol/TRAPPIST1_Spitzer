function histogram(param,nbin)
  p1 = minimum(param)-1e-15; p2 = maximum(param)+1e-15
  pbin_square = [p1]
  hist_square = [0.0]
  pbin = zeros(nbin)
  hist = zeros(nbin)
  psort = sort(param)
  i1 = 1; np = size(param)[1]
  println(p1," ",psort[i1]," ",p2," ",psort[np])
  dp = (p2-p1)/nbin
  for i=1:nbin
    while psort[i1] <= p1+dp*i
      hist[i] += 1.0
      i1 += 1
      if i1 == np+1
        break
      end
    end
    push!(hist_square,hist[i]); push!(hist_square,hist[i])
    pbin[i] = p1+dp*(i-0.5)
    push!(pbin_square,p1+dp*(i-1))
    push!(pbin_square,p1+dp*i)
    if i1 == np+1
      break
    end
  end
  push!(hist_square,0.0)
  push!(pbin_square,p2)
  return pbin,hist,pbin_square,hist_square
end
