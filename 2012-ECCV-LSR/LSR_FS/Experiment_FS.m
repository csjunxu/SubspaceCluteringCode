% Face clustering on the Extended Yale B database by using LSR in

% Can-Yi Lu, Hai Min, Zhong-Qiu Zhao, Lin Zhu, De-Shuang Huang and Shuicheng Yan.
% Robust and Efficient Subspace Segmentation via Least Squares Regression,
% European Conference on Computer Vision (ECCV), 2012.

%--------------------------------------------------------------------------
% Copyright @ Can-Yi Lu, 2012
%--------------------------------------------------------------------------


clear ;
close all;
currentpath = cd ;
AddedPath = genpath( currentpath ) ;
addpath( AddedPath ) ;
fprintf('\n\n**************************************   %s   *************************************\n' , datestr(now) );
fprintf( [ mfilename(currentpath) ' Begins.\n' ] ) ;
fprintf( [ mfilename(currentpath) ' is going, please wait...\n' ] ) ;


for nCluster = [5:1:10];           % number of subspace, 5 or 10 used in our paper
    %% Data YaleB
    load 'C:\Users\csjunxu\Desktop\SC\2012-ECCV-LSR\LSR_FS\Data\YaleB.mat'              % load YaleB dataset
    
    
    num = nCluster * 64 ;    % number of data used for subspace segmentation
    fea = fea(:,1:num) ;
    gnd = gnd(:,1:num) ;
    redDim = nCluster * 6 ;
    
    
    
    %% PCA Projection
    [ eigvector , eigvalue ] = PCA( fea ) ;
    maxDim = length(eigvalue)
    fea = eigvector' * fea ;
    
    % normalize
    for i = 1 : size(fea,2)
        fea(:,i) = fea(:,i) /norm(fea(:,i)) ;
    end
    
    
    %% Subspace segmentation methods
    % SegmentationMethod = 'LSR1' ;     % LSR1 by (16) in our paper
    SegmentationMethod = 'LSR2' ;     % LSR2 by (18) in our paper
    
    writefilepath = 'C:/Users/csjunxu/Desktop/SC/Results/';
    % writefilepath = '';
    
    %% Parameter
    
    switch nCluster
        case 5
            para = [0.4] * ones(1,20) ;
        case 10
            para = [0.004 ] * ones(1,20) ;
    end
    
    for lambda = [1e-6 5e-6 1e-5 5e-5 1e-4 5e-4 1e-3 5e-3 1e-2 5e-2 .1 .5 1 5]
        para = lambda * ones(1,20) ;
        %% Output results
        fid = 1 ;  % output to the screen
        fprintf( fid , ['Function                   = ' mfilename(currentpath) '.m\n'] ) ;
        fprintf( fid ,  'Data                       = %s, nCluster = %d\n' , 'YaleB' , nCluster ) ;
        fprintf( fid ,  'SegmentationMethod         = %s\n' , SegmentationMethod ) ;
        num_para = length( para ) ;
        num_redDim = length( redDim ) ;
        fprintf( fid , 'para =' ) ;
        for i = 1 : num_para
            fprintf( fid , '\t%5f' , para(i) ) ;
        end
        fprintf( fid , '\n' ) ;
        
        %% Subspace segmentation
        Accuracy = zeros( num_redDim , num_para ) ;
        for i = 1 : num_redDim
            d = redDim( i ) ;
            fprintf( 'd = %d', d ) ;
            Yfea = fea(1:d,:) ;
            for j = 1 : num_para
                p = para( j ) ;
                Accuracy(i,j) = SubspaceSegmentation( SegmentationMethod , Yfea , gnd , p ) ;
                fprintf( fid , '\t%.3f ' , Accuracy(i,j)*100 ) ;
            end
            fprintf('\n') ;
        end
        
        
        %% output
        fprintf('\n\n');
        fprintf( fid , 'para =' ) ;
        for i = 1 : length(para)
            fprintf( fid , '\t%5f' , para(i) ) ;
        end
        fprintf( fid , '\n' ) ;
        for i = 1 : length(redDim)
            d = redDim( i ) ;
            fprintf( 'd = %d', d ) ;
            for j = 1 : length(para)
                fprintf( fid , '\t%.3f ' , Accuracy(i,j)*100 ) ;
            end
            fprintf('\n') ;
        end
        
        [maxa ind] = max( Accuracy*100 )
        maxpara = para(ind)
        maxAcc = max( max(Accuracy*100) )
        
        %% output
        matname = sprintf([writefilepath 'YaleB_' SegmentationMethod '_DR' num2str(redDim) '_nCluster' num2str(nCluster) '_lambda' num2str(lambda) '.mat']);
        save(matname, 'Accuracy', 'maxAcc');
    end
end

% rmpath( genpath(currentpath) ) ;
