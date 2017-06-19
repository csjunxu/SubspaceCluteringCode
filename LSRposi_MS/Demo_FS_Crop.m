
clear ;

load 'C:\Users\csjunxu\Desktop\SC\Datasets\YaleB_Crop.mat'              % load YaleB dataset
writefilepath = 'C:/Users/csjunxu/Desktop/SC/Results/';
% writefilepath = '';

Repeat = 1; %number of repeations
DR = 1; % perform dimension reduction or not
dim = 6;

%% Subspace segmentation methods
% SegmentationMethod = 'LSR' ;
% SegmentationMethod = 'LSRd0' ;

% SegmentationMethod = 'NNLSR_LSR' ;
% SegmentationMethod = 'NNLSRd0_LSR' ;
% SegmentationMethod = 'NPLSR_LSR' ;
% SegmentationMethod = 'NPLSRd0_LSR' ;
% find a fast solver is still in process

% SegmentationMethod = 'ANNLSR_LSR' ;
% SegmentationMethod = 'ANNLSRd0_LSR' ;
SegmentationMethod = 'ANPLSR_LSR' ;
% SegmentationMethod = 'ANPLSRd0_LSR' ;

%% Subspace segmentation
for maxIter = [5 10]
    Par.maxIter = maxIter;
    for rho = [0.05 0.06 0.04]
        Par.rho = rho;
        for lambda = [4 5 6 7]
            Par.lambda = 10^(-lambda);
            for nSet = [2 3 5 8 10]
                n = nSet;
                index = Ind{n};
                for i = 1:size(index,1)
                    fea = [];
                    gnd = [];
                    for p = 1:n
                        fea = [fea Y{index(i, p), 1}];
                        gnd= [gnd p * ones(1, length(S{index(i, p)}))];
                    end
                    [D, N] = size(fea);
                    fprintf( '%d: %d\n', size(index, 1), i ) ;
                    redDim = size(fea, 1);
                    if DR == 1
                        %% PCA Projection
                        [ eigvector , eigvalue ] = PCA( fea ) ;
                        maxDim = length(eigvalue);
                        fea = eigvector' * fea ;
                        redDim = n * dim ;
                    end
                    %% normalize
                    for c = 1 : size(fea,2)
                        fea(:,c) = fea(:,c) /norm(fea(:,c)) ;
                    end
                    
                    %% Subspace Clustering
                    missrate = zeros(size(index, 1), Repeat) ;
                    fprintf( 'dimension = %d \n', redDim ) ;
                    Yfea = fea(1:redDim, :) ;
                    for j = 1 : Repeat
                        switch SegmentationMethod
                            case 'LSR'
                                C = LSR( Yfea , Par ) ;
                            case 'LSRd0'
                                C = LSRd0( Yfea , Par ) ; % solved by ADMM
                                % C = LSR1( ProjX , Par.lambda ) ; % proposed by Lu
                            case 'NNLSR_LSR'                   % non-negative
                                C = NNLSR( Yfea , Par ) ;
                            case 'NNLSRd0_LSR'               % non-negative, diagonal = 0
                                C = NNLSRd0( Yfea , Par ) ;
                            case 'NPLSR_LSR'                   % non-positive
                                C = NPLSR( Yfea , Par ) ;
                            case 'NPLSRd0_LSR'               % non-positive, diagonal = 0
                                C = NPLSRd0( Yfea , Par ) ;
                            case 'ANNLSR_LSR'                 % affine, non-negative
                                C = ANNLSR( Yfea , Par ) ;
                            case 'ANNLSRd0_LSR'             % affine, non-negative, diagonal = 0
                                C = ANNLSRd0( Yfea , Par ) ;
                            case 'ANPLSR_LSR'                 % affine, non-positive
                                C = ANPLSR( Yfea , Par ) ;
                            case 'ANPLSRd0_LSR'             % affine, non-positive, diagonal = 0
                                C = ANPLSRd0( Yfea , Par ) ;
                        end
                        for k = 1 : size(C,2)
                            C(:, k) = C(:, k) / max(abs(C(:, k))) ;
                        end
                        Z = ( abs(C) + abs(C') ) / 2 ;
                        idx = clu_ncut(Z,n) ;
                        missrate(i, j) = 1 - compacc(idx,gnd);
                        fprintf('%.3f%% \n' , missrate(i, j)*100) ;
                    end
                    missrateTot{n}(i) = mean(missrate(i, :)*100);
                    fprintf('Mean error of %d/%d is %.3f%%.\n ' , i, size(index, 1), missrateTot{n}(i)) ;
                end
                %% output
                avgmissrate(n) = mean(missrateTot{n});
                medmissrate(n) = median(missrateTot{n});
                fprintf('Total mean error  is %.3f%%.\n ' , avgmissrate(n)) ;
                matname = sprintf([writefilepath 'YaleB_Crop_' SegmentationMethod '_dim' num2str(dim) '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_lambda' num2str(Par.lambda) '.mat']);
                save(matname,'missrateTot','avgmissrate','medmissrate');
            end
        end
    end
end
