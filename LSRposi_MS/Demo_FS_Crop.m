
clear ;

load 'C:\Users\csjunxu\Desktop\SC\Datasets\YaleB_Crop.mat'              % load YaleB dataset
writefilepath = 'C:/Users/csjunxu/Desktop/SC/Results/';
% writefilepath = '';



Repeat = 20; %number of repeations
DR = 1; % perform dimension reduction or not


%% Data YaleB
nSet = [2:1:10];
for set = 1:length(nSet)
    n = nSet(set);
    index = Ind{n};
    for j = 1:size(index,1)
        fea = [];
        gnd = [];
        for p = 1:n
            fea = [fea Y{p, 1}];
            gnd= [gnd p * ones(1, length(S{p}))];
        end
        [D, N] = size(fea);
        
        redDim = size(fea, 1);
        if DR == 1
            %% PCA Projection
            [ eigvector , eigvalue ] = PCA( fea ) ;
            maxDim = length(eigvalue);
            fea = eigvector' * fea ;
            redDim = nSet * 6 ;
        end
        %% normalize
        for i = 1 : size(fea,2)
            fea(:,i) = fea(:,i) /norm(fea(:,i)) ;
        end
        
        %% Subspace segmentation methods
        % SegmentationMethod = 'LSRd0po_LSR' ;
        SegmentationMethod = 'LSRpo_LSR' ;
        % SegmentationMethod = 'LSRd0ne_LSR' ;
        % SegmentationMethod = 'LSRne_LSR' ;
        % SegmentationMethod = 'LSRd0_LSR' ;
        
        %% Output results
        num_redDim = length( redDim ) ;
        
        %% Subspace segmentation
        for maxIter = [2 5 10 15 20 25 30]
            Par.maxIter = maxIter;
            for lambda = [.00001 .00005 .0001 .0005 .001:.001:.01]
                Par.lambda = lambda;
                for rho = [0.02]
                    Par.rho = rho;
                    Accuracy = zeros( num_redDim , Repeat ) ;
                    for i = 1 : num_redDim
                        d = redDim( i ) ;
                        fprintf( 'd = %d', d ) ;
                        Yfea = fea(1:d, :) ;
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
                            nSet = length( unique( gnd ) ) ;
                            Z = ( abs(C) + abs(C') ) / 2 ;
                            idx = clu_ncut(Z,nSet) ;
                            Accuracy(i,j) = compacc(idx,gnd);
                            fprintf( fid , '\t%.3f ' , Accuracy(i,j)*100 ) ;
                        end
                        fprintf('\n') ;
                    end
                    %% output
                    fprintf('\n\n');
                    for i = 1 : length(redDim)
                        d = redDim( i ) ;
                        fprintf( 'd = %d', d ) ;
                        for j = 1 : Repeat
                            fprintf( '%.3f ' , Accuracy(i,j)*100 ) ;
                        end
                        fprintf('\n') ;
                    end
                    [maxa ind] = max( Accuracy*100 )
                    maxAcc = max( max(Accuracy*100) )
                    %% output
                    matname = sprintf([writefilepath 'YaleB_Crop_' SegmentationMethod '_DR' num2str(redDim) '_nCluster' num2str(nSet) '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_lambda' num2str(lambda) '.mat']);
                    save(matname, 'Accuracy', 'maxAcc');
                end
            end
        end
    end
end
