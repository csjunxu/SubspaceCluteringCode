
clear ;

load 'C:\Users\csjunxu\Desktop\SC\Datasets\YaleBCrop025.mat';
% load 'C:\Users\csjunxu\Desktop\SC\Datasets\USPS_Crop.mat'   % load USPS dataset


writefilepath = 'C:/Users/csjunxu/Desktop/SC/Results/';

%% Subspace segmentation methods
% SegmentationMethod = 'LSRd0po_SSC' ;
SegmentationMethod = 'LSRpo_SSC' ;
% SegmentationMethod = 'LSRd0ne_SSC' ;
% SegmentationMethod = 'LSRne_SSC' ;
% SegmentationMethod = 'LSRd0_SSC' ;
DR = 1; % dimension reduction
dim = 6;
%% Subspace segmentation
for maxIter = [10 15]
    Par.maxIter = maxIter;
    for mu = [1]
        Par.mu = mu;
        for rho = [0.1:0.1:0.5]
            Par.rho = rho;
            for lambda = [2:1:6]
                Par.lambda = 10^(-lambda);
                for nSet = [2 3 5 8 10]
                    for i = 1:length(nSet)
                        n = nSet(i);
                        index = Ind{n};
                        for j = 1:size(index,1)
                            X = [];
                            for p = 1:n
                                X = [X Y(:,:,index(j,p))];
                            end
                            [D,N] = size(X);
                            
                            fea = X ;
                            gnd = s{n} ;
                            
                            redDim = size(fea, 1);
                            if DR == 1
                                %% PCA Projection
                                [ eigvector , eigvalue ] = PCA( fea ) ;
                                maxDim = length(eigvalue);
                                fea = eigvector' * fea ;
                                redDim = nSet * dim ;
                            end
                            
                            % normalize
                            for c = 1 : size(fea,2)
                                fea(:,c) = fea(:,c) /norm(fea(:,c)) ;
                            end
                            Yfea = fea(1:redDim, :) ;
                            switch SegmentationMethod
                                case 'LSRd0po_SSC'
                                    C = LSRd0po( Yfea , Par ) ;
                                case 'LSRpo_SSC'
                                    C = LSRpo( Yfea , Par ) ;
                                case 'LSRd0ne_SSC'
                                    C = LSRd0ne( Yfea , Par ) ;
                                case 'LSRne_SSC'
                                    C = LSRne( Yfea , Par ) ;
                                case 'LSRd0_SSC'
                                    C = LSRd0( Yfea , Par ) ;
                            end
                            for k = 1 : size(C,2)
                                C(:, k) = C(:, k) / max(abs(C(:, k))) ;
                            end
                            nCluster = length( unique( gnd ) ) ;
                            Z = ( abs(C) + abs(C') ) / 2 ;
                            idx = clu_ncut(Z,nCluster) ;
                            accuracy = compacc(idx,gnd');
                            missrate = 1 - accuracy;
                            fprintf('missrate: %2.4f\n', missrate);
                            missrateTot{n}(j) = missrate;
                        end
                        avgmissrate(n) = mean(missrateTot{n});
                        medmissrate(n) = median(missrateTot{n});
                        matname = sprintf([writefilepath 'YaleB_' SegmentationMethod '_DR' num2str(DR) '_dim' num2str(dim) '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_mu' num2str(Par.mu) '_lambda' num2str(Par.lambda) '.mat']);
                        save(matname,'missrateTot','avgmissrate','medmissrate');
                    end
                    matname = sprintf([writefilepath 'YaleB_' SegmentationMethod '_DR' num2str(DR) '_dim' num2str(dim) '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_mu' num2str(Par.mu) '_lambda' num2str(Par.lambda) '.mat']);
                    save(matname,'missrateTot','avgmissrate','medmissrate');
                end
            end
        end
    end
end




