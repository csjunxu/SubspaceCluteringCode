
clear ;
close all;
currentpath = cd ;
AddedPath = genpath( currentpath ) ;
addpath( AddedPath ) ;
fprintf('\n\n**************************************   %s   *************************************\n' , datestr(now) );
fprintf( [ mfilename(currentpath) ' Begins.\n' ] ) ;
fprintf( [ mfilename(currentpath) ' is going, please wait...\n' ] ) ;

%% Data YaleB
load 'C:\Users\csjunxu\Desktop\SC\2012-ECCV-LSR\LSR_FS\Data\YaleB.mat'              % load YaleB dataset
nCluster = 5 ;           % number of subspace, 5 or 10 used in our paper
num = nCluster * 64 ;    % number of data used for subspace segmentation
fea = fea(:,1:num) ;
gnd = gnd(:,1:num) ;

% writefilepath = 'C:/Users/csjunxu/Desktop/SC/Results/';
writefilepath = '';
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

%% Output results
fid = 1 ;  % output to the screen
fprintf( fid , ['Function                   = ' mfilename(currentpath) '.m\n'] ) ;
fprintf( fid ,  'Data                       = %s, nCluster = %d\n' , 'YaleB' , nCluster ) ;
fprintf( fid ,  'SegmentationMethod         = %s\n' , SegmentationMethod ) ;
fprintf( fid , '\n' ) ;

%% Subspace segmentation
for maxIter = [10 15 20 25 30]
    Par.maxIter = maxIter;
    for mu = [1]
        Par.mu = mu;
        for lambda = [.1 .11 .12 .09 .08 .13 .14 .07 .06 .15 .05]
            Par.lambda = lambda;
            for rho = [0.02]
                Par.rho = rho;
                    Yfea = fea(1:redDim,:) ;
                switch SegmentationMethod
                    case 'LSRd0po_LSR'
                        C = LSRd0po( Yfea , Par ) ;
                    case 'LSRpo_LSR'
                        C = LSRpo( Yfea , Par ) ;
                end
                
                for i = 1 : size(C,2)
                    C(:,i) = C(:,i) / max(abs(C(:,i))) ;
                end
                nCluster = length( unique( gnd ) ) ;
                Z = ( abs(C) + abs(C') ) / 2 ;
                idx = clu_ncut(Z,nCluster) ;
                Accuracy = compacc(idx,gnd);
                fprintf( fid , '\t%.3f \% \n' , Accuracy*100 ) ;
                
                %% output
                matname = sprintf([writefilepath 'YaleB_' SegmentationMethod '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_mu' num2str(Par.mu) '_lambda' num2str(lambda) '.mat']);
                save(matname,'Accuracy');
            end
        end
    end
end



currentpath = cd ;
AddedPath = genpath( currentpath ) ;
addpath( AddedPath ) ;
fprintf('\n\n**************************************   %s   *************************************\n' , datestr(now) );
fprintf( [ mfilename(currentpath) ' Begins.\n' ] ) ;
fprintf( [ mfilename(currentpath) ' is going, please wait...\n' ] ) ;

%% Data YaleB
load 'C:\Users\csjunxu\Desktop\SC\2012-ECCV-LSR\LSR_FS\Data\YaleB.mat'              % load YaleB dataset
nCluster = 10 ;           % number of subspace, 5 or 10 used in our paper
num = nCluster * 64 ;    % number of data used for subspace segmentation
fea = fea(:,1:num) ;
gnd = gnd(:,1:num) ;
redDim = nCluster * 6 ;

%% PCA Projection
[ eigvector , eigvalue ] = PCA( fea ) ;
maxDim = length(eigvalue);
fea = eigvector' * fea ;

% normalize
for i = 1 : size(fea,2)
    fea(:,i) = fea(:,i) /norm(fea(:,i)) ;
end

%% Output results
fid = 1 ;  % output to the screen
fprintf( fid , ['Function                   = ' mfilename(currentpath) '.m\n'] ) ;
fprintf( fid ,  'Data                       = %s, nCluster = %d\n' , 'YaleB' , nCluster ) ;
fprintf( fid ,  'SegmentationMethod         = %s\n' , SegmentationMethod ) ;
fprintf( fid , '\n' ) ;

%% Subspace segmentation
for maxIter = [10 15 20 25 30]
    Par.maxIter = maxIter;
    for mu = [1]
        Par.mu = mu;
        for lambda = [.1 .11 .12 .09 .08 .13 .14 .07 .06 .15 .05]
            Par.lambda = lambda;
            for rho = [0.02]
                Par.rho = rho;
            Yfea = fea(1:redDim,:) ;
                
                switch SegmentationMethod
                    case 'LSRd0po_LSR'
                        C = LSRd0po( Yfea , Par ) ;
                    case 'LSRpo_LSR'
                        C = LSRpo( Yfea , Par ) ;
                end
                
                for i = 1 : size(C,2)
                    C(:,i) = C(:,i) / max(abs(C(:,i))) ;
                end
                nCluster = length( unique( gnd ) ) ;
                Z = ( abs(C) + abs(C') ) / 2 ;
                idx = clu_ncut(Z,nCluster) ;
                Accuracy = compacc(idx,gnd);
                fprintf( fid , '\t%.3f \% \n' , Accuracy*100 ) ;
                
                %% output
                matname = sprintf([writefilepath 'YaleB_' SegmentationMethod '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_mu' num2str(Par.mu) '_lambda' num2str(lambda) '.mat']);
                save(matname,'Accuracy');
            end
        end
    end
end



