%----------------------------------------------------------------------------------------------------
% Copyright @ Jun Xu, 2017
% Please Contact csjunxu@comp.polyu.edu.hk for any questions.
%----------------------------------------------------------------------------------------------------

clc, clear all, close all
addpath(genpath(fullfile(pwd)));

% Please download Hopkins155 dataset at http://www.vision.jhu.edu/data/hopkins155/

cd 'C:/Users/csjunxu/Desktop/SC/Datasets/Hopkins155/';

method = 'LSRd0po_SSC';
% method = 'LSRpo_SSC';
for mu = [1]
    Par.mu = mu;
    for maxIter = [200]
        Par.maxIter = maxIter;
        for lambda = [2e-4:1e-4:1e-3]
            Par.lambda = lambda;
            for rho = [0.1 0.01 0.001]
                Par.rho = rho;
                
                maxNumGroup = 5;
                for i = 1:maxNumGroup
                    num(i) = 0;
                end
                
                d = dir;
                for i = 1:length(d)
                    if ( (d(i).isdir == 1) && ~strcmp(d(i).name,'.') && ~strcmp(d(i).name,'..') )
                        filepath = d(i).name;
                        eval(['cd ' filepath]);
                        
                        f = dir;
                        foundValidData = false;
                        for j = 1:length(f)
                            if ( ~isempty(strfind(f(j).name,'_truth.mat')) )
                                ind = j;
                                foundValidData = true;
                                break
                            end
                        end
                        eval(['load ' f(ind).name]);
                        cd ..
                        
                        if (foundValidData)
                            n = max(s);
                            N = size(x,2);
                            F = size(x,3);
                            D = 2*F;
                            X = reshape(permute(x(1:2,:,:),[1 3 2]),D,N);
                            
                            switch SegmentationMethod
                                case 'LSRd0po_SSC'
                                    r = 0; affine = false; outlier = false; rho = 0.7;
                                    missrate1= LSRd0posi(X,r,affine,outlier,rho,s,Par);
                                    
                                    r = 4*n; affine = false; outlier = false; rho = 0.7;
                                    missrate2 = LSRd0posi(X,r,affine,outlier,rho,s,Par);
                                case 'LSRpo_SSC'
                                    r = 0; affine = false; outlier = false; rho = 0.7;
                                    missrate1= LSRposi(X,r,affine,outlier,rho,s,Par);
                                    
                                    r = 4*n; affine = false; outlier = false; rho = 0.7;
                                    missrate2 = LSRposi(X,r,affine,outlier,rho,s,Par);
                            end
                            
                            fprintf('%d %s\t %f, %f\n', i, f(ind).name , missrate1, missrate2 ) ;
                            num(n) = num(n) + 1;
                            missrateTot1{n}(num(n)) = missrate1;
                            missrateTot2{n}(num(n)) = missrate2;
                            
                            eval(['cd ' filepath]);
                            cd ..
                        end
                    end
                end
                
                L = [2 3];
                allmissrate1 = [];
                allmissrate2 = [];
                for i = 1:length(L)
                    j = L(i);
                    avgmissrate1(j) = mean(missrateTot1{j});
                    medmissrate1(j) = median(missrateTot1{j});
                    allmissrate1 = [allmissrate1 missrateTot1{j}];
                    
                    avgmissrate2(j) = mean(missrateTot2{j});
                    medmissrate2(j) = median(missrateTot2{j});
                    allmissrate2 = [allmissrate2 missrateTot2{j}];
                end
                avgallmissrate1 = sum(allmissrate1)/length(allmissrate1);
                medallmissrate1 = median(allmissrate1);
                avgallmissrate2 = sum(allmissrate2)/length(allmissrate2);
                medallmissrate2 = median(allmissrate2);
                matname = sprintf(['C:/Users/csjunxu/Desktop/SC/Results/' method '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_mu' num2str(Par.mu) '_lambda' num2str(lambda) '.mat']);
                save(matname,'avgallmissrate1','medallmissrate1','missrateTot1','avgmissrate1','medmissrate1','avgallmissrate2','medallmissrate2','missrateTot2','avgmissrate2','medmissrate2');
            end
        end
    end
end