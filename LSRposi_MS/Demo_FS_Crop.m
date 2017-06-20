
clear ;

load 'C:\Users\csjunxu\Desktop\SC\Datasets\YaleB_Crop.mat'              % load YaleB dataset
writefilepath = 'C:/Users/csjunxu/Desktop/SC/Results/';
% writefilepath = '';
dataset = 'YaleB_LSR';
Repeat = 1; %number of repeations
DR = 1; % perform dimension reduction or not
dim = 6;

%% Subspace segmentation methods
% SegmentationMethod = 'LSR' ;
% SegmentationMethod = 'LSRd0' ;
SegmentationMethod = 'LSR1' ;
% SegmentationMethod = 'LSR2' ;

% SegmentationMethod = 'NNLSR' ;
% SegmentationMethod = 'NNLSRd0' ;
% SegmentationMethod = 'NPLSR' ;
% SegmentationMethod = 'NPLSRd0' ;
% find a fast solver is still in process

% SegmentationMethod = 'ANNLSR' ;
% SegmentationMethod = 'ANNLSRd0' ;
% SegmentationMethod = 'ANPLSR' ;
% SegmentationMethod = 'ANPLSRd0' ;

%% Subspace segmentation
for maxIter = [5]
    Par.maxIter = maxIter;
    for rho = [0.05]
        Par.rho = rho;
        for lambda = [1:1:5]
            Par.lambda = lambda*10^(-3);
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
                            case 'LSR1'
                                C = LSR1( Yfea , Par.lambda ) ; % proposed by Lu
                            case 'LSR2'
                                C = LSR2( Yfea , Par.lambda ) ; % proposed by Lu
                            case 'LSR'
                                C = LSR( Yfea , Par ) ;
                            case 'LSRd0'
                                C = LSRd0( Yfea , Par ) ; % solved by ADMM
                            case 'NNLSR'                   % non-negative
                                C = NNLSR( Yfea , Par ) ;
                            case 'NNLSRd0'               % non-negative, diagonal = 0
                                C = NNLSRd0( Yfea , Par ) ;
                            case 'NPLSR'                   % non-positive
                                C = NPLSR( Yfea , Par ) ;
                            case 'NPLSRd0'               % non-positive, diagonal = 0
                                C = NPLSRd0( Yfea , Par ) ;
                            case 'ANNLSR'                 % affine, non-negative
                                C = ANNLSR( Yfea , Par ) ;
                            case 'ANNLSRd0'             % affine, non-negative, diagonal = 0
                                C = ANNLSRd0( Yfea , Par ) ;
                            case 'ANPLSR'                 % affine, non-positive
                                C = ANPLSR( Yfea , Par ) ;
                            case 'ANPLSRd0'             % affine, non-positive, diagonal = 0
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
                matname = sprintf([writefilepath dataset '_' SegmentationMethod '_DR' num2str(DR) '_dim' num2str(dim) '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_lambda' num2str(Par.lambda) '.mat']);
                save(matname,'missrateTot','avgmissrate','medmissrate');
            end
        end
    end
end
