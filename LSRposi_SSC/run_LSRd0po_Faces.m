clear all;

load 'C:\Users\csjunxu\Desktop\SC\Datasets\YaleBCrop025.mat';

method = 'LSRd0po_SSC';
% method = 'LSRpo_SSC';

for mu = [1]
    Par.mu = mu;
    for maxIter = [150 100 200]
        Par.maxIter = maxIter;
        for lambda = [0.1 0.5 1]
            Par.lambda = lambda;
            for rho = [0.01 0.1]
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
                        matname = sprintf(['C:/Users/csjunxu/Desktop/SC/Results/YaleB_' method '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_mu' num2str(Par.mu) '_lambda' num2str(lambda) '.mat']);
                        save(matname,'avgallmissrate1','medallmissrate1','missrateTot1','avgmissrate1','medmissrate1','avgallmissrate2','medallmissrate2','missrateTot2','avgmissrate2','medmissrate2');
                    end
                    matname = sprintf(['C:/Users/csjunxu/Desktop/SC/Results/YaleB_' method '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_mu' num2str(Par.mu) '_lambda' num2str(lambda) '.mat']);
                    save(matname,'avgallmissrate1','medallmissrate1','missrateTot1','avgmissrate1','medmissrate1','avgallmissrate2','medallmissrate2','missrateTot2','avgmissrate2','medmissrate2');
                end
            end
        end
    end
end