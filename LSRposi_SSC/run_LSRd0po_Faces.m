% error is too high
% possible solutions
% 1. originally, the number of iteration may not be enough at 200.
% 2. Simulate the settings of SSC

clear all;

load 'C:\Users\csjunxu\Desktop\SC\Datasets\YaleBCrop025.mat';

method = 'LSRd0po_SSC';
% method = 'LSRpo_SSC';
Par.thr = 2*10^-4;
for mu = [1]
    Par.mu = mu;
    for maxIter = 200
        Par.maxIter = maxIter;
        for rho = [0.5]
            Par.rho = rho;
            for nSet = 2 %[2 3 4 5 6 7 8 9 10]
                for i = 1:length(nSet)
                    n = nSet(i);
                    idx = Ind{n};
                    for j = 1:size(idx,1)
                        X = [];
                        for p = 1:n
                            X = [X Y(:,:,idx(j,p))];
                        end
                        [D,N] = size(X);
                        
                        r = 0; affine = false; outlier = true; rho = 1;
                        time0  =   clock;
                        missrate= LSRd0posi(X,r,affine,outlier,rho,s{n},Par);
                        missrateTot{n}(j) = missrate;
                    end
                    avgmissrate(n) = mean(missrateTot{n});
                    medmissrate(n) = median(missrateTot{n});
                    matname = sprintf(['C:/Users/csjunxu/Desktop/SC/Results/YaleB_' method '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_mu' num2str(Par.mu) '.mat']);
                    save(matname,'missrateTot','avgmissrate','medmissrate');
                end
                matname = sprintf(['C:/Users/csjunxu/Desktop/SC/Results/YaleB_' method '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_mu' num2str(Par.mu) '.mat']);
                save(matname,'missrateTot','avgmissrate','medmissrate');
            end
        end
    end
end