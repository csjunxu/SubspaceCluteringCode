
clear ;
close all;
currentpath = cd ;
AddedPath = genpath( currentpath ) ;
addpath( AddedPath ) ;
fprintf('\n\n**************************************   %s   *************************************\n' , datestr(now) );
fprintf( [ mfilename(currentpath) ' Begins.\n' ] ) ;
fprintf( [ mfilename(currentpath) ' is going, please wait...\n' ] ) ;

%% Data YaleB
% load 'C:\Users\csjunxu\Desktop\SC\Datasets\YaleBCrop025.mat';
load 'C:\Users\csjunxu\Desktop\SC\2012-ECCV-LSR\LSR_FS\Data\YaleB.mat'              % load YaleB dataset
nCluster = 10 ;           % number of subspace, 5 or 10 used in our paper
num = nCluster * 64 ;    % number of data used for subspace segmentation
fea = fea(:,1:num) ;
gnd = gnd(:,1:num) ;

writefilepath = 'C:/Users/csjunxu/Desktop/SC/Results/';
% writefilepath = '';
%% PCA Projection
[ eigvector , eigvalue ] = PCA( fea ) ;
maxDim = length(eigvalue);
fea = eigvector' * fea ;
redDim = nCluster * 6 ;

% normalize
for i = 1 : size(fea,2)
    fea(:,i) = fea(:,i) /norm(fea(:,i)) ;
end

%% Subspace segmentation methods
% SegmentationMethod = 'LSRd0po_LSR' ;
SegmentationMethod = 'LSRpo_LSR' ;
% SegmentationMethod = 'LSRd0ne_LSR' ;
% SegmentationMethod = 'LSRne_LSR' ;
% SegmentationMethod = 'LSRd0_LSR' ;

%% Parameter
switch nCluster
    case 5
        Repeat = 20;
    case 10
        Repeat = 20;
end


%% Output results
fid = 1 ;  % output to the screen
fprintf( fid , ['Function                   = ' mfilename(currentpath) '.m\n'] ) ;
fprintf( fid ,  'Data                       = %s, nCluster = %d\n' , 'YaleB' , nCluster ) ;
fprintf( fid ,  'SegmentationMethod         = %s\n' , SegmentationMethod ) ;
num_redDim = length( redDim ) ;
fprintf( fid , '\n' ) ;

%% Subspace segmentation
for maxIter = [5 10 15 20 25 30]
    Par.maxIter = maxIter;
    for mu = [1]
        Par.mu = mu;
        for lambda = [.0009:-.0001:.0001]
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
                        nCluster = length( unique( gnd ) ) ;
                        Z = ( abs(C) + abs(C') ) / 2 ;
                        idx = clu_ncut(Z,nCluster) ;
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
                        fprintf( fid , '\t%.3f ' , Accuracy(i,j)*100 ) ;
                    end
                    fprintf('\n') ;
                end
                [maxa ind] = max( Accuracy*100 )
                maxAcc = max( max(Accuracy*100) )
                %% output
                matname = sprintf([writefilepath 'YaleB_LSR_' SegmentationMethod '_nCluster' num2str(nCluster) '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_mu' num2str(Par.mu) '_lambda' num2str(lambda) '.mat']);
                save(matname, 'Accuracy', 'maxAcc');
            end
        end
    end
end

