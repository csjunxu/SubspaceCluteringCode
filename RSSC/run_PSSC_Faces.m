clear;

load 'C:\Users\csjunxu\Desktop\SC\Datasets\YaleBCrop025.mat';

for alpha = 20;
    for nSet = [5]% [2 3 4 5 6 7 8 9 10]
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
                missrate = PSSC(X,r,affine,alpha,outlier,rho,s{n});
                fprintf('missrate: %2.4f\n', missrate);
                missrateTot{n}(j) = missrate;
            end
            avgmissrate(n) = mean(missrateTot{n});
            medmissrate(n) = median(missrateTot{n});
            save(['PSSC_Faces' num2str(alpha) '.mat'], missrateTot, avgmissrate, medmissrate, alpha);
        end
        save(['PSSC_Faces' num2str(alpha) '.mat'], missrateTot, avgmissrate, medmissrate, alpha);
    end
end