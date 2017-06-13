function missrate = LSRd0posi(X,r,affine,outlier,rho,s,Par)

if (nargin < 7)
    Par.lambda = 1;
    Par.rho = 0.1;
    Par.maxIter = 200;
end
if (nargin < 5)
    rho = 1;
end
if (nargin < 4)
    outlier = false;
end
if (nargin < 3)
    affine = false;
end
if (nargin < 2)
    r = 0;
end
n = max(s);
if r==0
    Xp = X ;
else
    [ eigvector , ~ ] = PCA( X ) ;
    Xp = eigvector(:,1:r)' * X ;
end

if (~outlier)
    CMat = LSRd0posi_admm(Xp,affine,Par);
    C = CMat;
else
%     CMat = LSRd0posi_admmOutlier(Xp,affine,Par);
CMat = admmOutlier_mat_func(Xp,affine,Par);
    N = size(Xp,2);
    C = CMat(1:N,:);
end
CKSym = BuildAdjacency(thrC(C,rho));
grps = SpectralClustering(CKSym,n);
missrate = compacc_ce(grps,s);