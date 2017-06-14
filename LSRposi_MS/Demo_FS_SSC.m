
clear ;

% error is too high
% possible solutions
% 1. originally, the number of iteration may not be enough at 200.
% 2. Simulate the settings of SSC

clear all;

load 'C:\Users\csjunxu\Desktop\SC\Datasets\YaleBCrop025.mat';

writefilepath = 'C:/Users/csjunxu/Desktop/SC/Results/';



%% Subspace segmentation methods
% SegmentationMethod = 'LSRd0po_LSR' ;
SegmentationMethod = 'LSRpo_LSR' ;
% SegmentationMethod = 'LSRd0ne_LSR' ;
% SegmentationMethod = 'LSRne_LSR' ;
% SegmentationMethod = 'LSRd0_LSR' ;
%% Subspace segmentation
for maxIter = [5 10 15 20 25 30]
    Par.maxIter = maxIter;
    for mu = [1]
        Par.mu = mu;
        for lambda = [.0009:-.0001:.0001]
            Par.lambda = lambda;
            for rho = [0.02]
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

                            % load 'C:\Users\csjunxu\Desktop\SC\2012-ECCV-LSR\LSR_FS\Data\YaleB.mat'              % load YaleB dataset
                            % nCluster = 10 ;           % number of subspace, 5 or 10 used in our paper
                            % num = nCluster * 64 ;    % number of data used for subspace segmentation
                            % fea = fea(:,1:num) ;
                            % gnd = gnd(:,1:num) ;

                            %                             %% PCA Projection
                            %                             [ eigvector , eigvalue ] = PCA( fea ) ;
                            %                             maxDim = length(eigvalue);
                            %                             fea = eigvector' * fea ;
                            %                             redDim = nCluster * 6 ;
                            
                            % normalize
                            for i = 1 : size(fea,2)
                                fea(:,i) = fea(:,i) /norm(fea(:,i)) ;
                            end
                            Yfea = fea;
                            switch SegmentationMethod
                                case 'LSRd0po_LSR'
                                    C = LSRd0po( Yfea , Par ) ;
                                case 'LSRpo_LSR'
                                    C = LSRpo( Yfea , Par ) ;
                                case 'LSRd0ne_LSR'
                                    C = LSRd0ne( Yfea , Par ) ;
                                case 'LSRne_LSR'
                                    C = LSRne( Yfea , Par ) ;
                                case 'LSRd0_LSR'
                                    C = LSRd0( Yfea , Par ) ;
                            end
                            for k = 1 : size(C,2)
                                C(:, k) = C(:, k) / max(abs(C(:, k))) ;
                            end
                            nCluster = length( unique( gnd ) ) ;
                            Z = ( abs(C) + abs(C') ) / 2 ;
                            idx = clu_ncut(Z,nCluster) ;
                            Accuracy = compacc(idx,gnd);
                            fprintf( fid , '\t%.3f ' , Accuracy*100 ) ;
                        end
                        fprintf('\n') ;
                        
                        fprintf('missrate: %2.4f\n', missrate);
                        missrateTot{n}(j) = missrate;
                    end
                    avgmissrate(n) = mean(missrateTot{n});
                    medmissrate(n) = median(missrateTot{n});
                    matname = sprintf([writefilepath 'YaleB_' method '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_mu' num2str(Par.mu) '.mat']);
                    save(matname,'missrateTot','avgmissrate','medmissrate');
                end
                matname = sprintf([writefilepath 'YaleB_' method '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_mu' num2str(Par.mu) '.mat']);
                save(matname,'missrateTot','avgmissrate','medmissrate');
            end
        end
    end
end




