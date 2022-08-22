function score = MAGComp_MVGML_Score(MAG,idxParents,CovMat,nSamples,tol)
% ML score single MAG district / c-component using RICF (based on GSMAG,Tria2016)
% note: slightly confusing use of 'component' and 'district'
score = 0;
% 
N = size(MAG,1);
if (N == 1) 
  % single node component?
  score = logdet(2*pi*CovMat)+(nSamples-1)/nSamples;    
  % note: GSMAG uses (n-1)/n where Nowzohour has n/(n+1) ...
else
  % multi-node component: fit and score
  % convert our MAG representation to other (transpose + tail mark == 3
  compMAG = MAG'; compMAG(compMAG == 1) = 3;
  % now call the RICF-fit routine, where we are only interested in the
  % implied (optimal) covariance matrix 'HatCovMat' given the MAG structure
  [~, ~, HatCovMat, ~] = GPS_RICF_fit(compMAG, CovMat, tol);
  % NOTE: RICF can often fail to converge for larger districts .. check 
  
  % get nr. of district nodes in component (=everything except the parents)
  Parents = find(idxParents == 1);  % parents in component MAG
  nrPa = length(Parents);
  compSize = N - nrPa;  % = |C_k|

  if any(idxParents)
    l1 = compSize*log(2*pi);
    l2 = logdet(HatCovMat) - log(prod(diag(HatCovMat(Parents, Parents))));
    l3 = (nSamples-1)/nSamples*(trace(HatCovMat\CovMat)-nrPa);
    score = l1+l2+l3;
  else
    l1 = compSize*log(2*pi);
    l2 = logdet(HatCovMat);
    l3 = (nSamples-1)/nSamples*trace(HatCovMat\CovMat);
    score = l1+l2+l3;
  end
    
end;  % if

end  % function Comp_MVGML_Score