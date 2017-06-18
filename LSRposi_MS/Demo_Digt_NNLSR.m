
clear ;

% load 'C:\Users\csjunxu\Desktop\SC\Datasets\YaleB_Crop.mat'  % load YaleB dataset
load 'C:\Users\csjunxu\Desktop\SC\Datasets\USPS_Crop.mat'   % load USPS dataset
% load 'C:\Users\csjunxu\Desktop\SC\Datasets\MNIST_Crop.mat' % load MNIST dataset


writefilepath = 'C:/Users/csjunxu/Desktop/SC/Results/';
% writefilepath = '';

Repeat = 1; %number of repeations
DR = 1; % perform dimension reduction or not
dim = 6;

%% Subspace segmentation
for maxIter = [5 10 15 20]
    Par.maxIter = maxIter;
    for rho = [0.05 0.06 0.04 0.07 0.03]
        Par.rho = rho;
        for lambda = [2:1:6]
            Par.lambda = 10^(-lambda);
            for nSet = [2 3 5 8 10] % nSet = [2:1:10];
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
                    
                    %% Subspace segmentation methods
                    % SegmentationMethod = 'LSRd0po_LSR' ;
                    SegmentationMethod = 'LSRpo_LSR' ;
                    % SegmentationMethod = 'LSRd0ne_LSR' ;
                    % SegmentationMethod = 'LSRne_LSR' ;
                    % SegmentationMethod = 'LSRd0_LSR' ;
                    
                    %% Subspace Clustering
                    missrate = zeros(size(index, 1), Repeat) ;
                    fprintf( 'dimension = %d \n', redDim ) ;
                    Yfea = fea(1:redDim, :) ;
                    for j = 1 : Repeat
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
                        Z = ( abs(C) + abs(C') ) / 2 ;
                        idx = clu_ncut(Z,n) ;
                        missrate(i, j) = 1 - compacc(idx,gnd);
                        fprintf('%.3f%% \n' , missrate(i, j)*100) ;
                    end
                    missrateTot{n}(i) = mean(missrate(i, :)*100);
                    fprintf('Mean Accuracy of %d/%d is %.3f%%.\n ' , i, size(index, 1), missrateTot{n}(i)) ;
                end
                %% output
                avgmissrate(n) = mean(missrateTot{n});
                medmissrate(n) = median(missrateTot{n});
                fprintf('Total mean missrate  is %.3f%%.\n ' , avgmissrate(n)) ;
                matname = sprintf([writefilepath 'USPS_Crop_' SegmentationMethod '_DR' num2str(redDim) '_dim' num2str(dim) '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_lambda' num2str(Par.lambda) '.mat']);
                save(matname,'missrateTot','avgmissrate','medmissrate');
            end
        end
    end
end
