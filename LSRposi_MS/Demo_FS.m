
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
redDim = nCluster * 6 ;



%% PCA Projection
[ eigvector , eigvalue ] = PCA( fea ) ;
maxDim = length(eigvalue);
fea = eigvector' * fea ;

% normalize
for i = 1 : size(fea,2)
    fea(:,i) = fea(:,i) /norm(fea(:,i)) ;
end


%% Subspace segmentation methods
SegmentationMethod = 'LSRd0po' ;
% SegmentationMethod = 'LSRpo' ;



%% Parameter

switch nCluster
    case 5
        para = [0.4] * ones(1,20) ;
    case 10
        para = [0.004 ] * ones(1,20) ;
end



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

for mu = [1]
    Par.mu = mu;
    for lambda = [0.1 0.5 1]
        Par.lambda = lambda;
        for rho = [0.1 0.01 0.05]
            Par.rho = rho;
            for maxIter = [200]
                Par.maxIter = maxIter;
                
                Accuracy = zeros( num_redDim , num_para ) ;
                for i = 1 : num_redDim
                    d = redDim( i ) ;
                    fprintf( 'd = %d', d ) ;
                    Yfea = fea(1:d,:) ;
                    for j = 1 : num_para
                        Accuracy(i,j) = SubspaceSegmentation( SegmentationMethod , Yfea , gnd , Par ) ;
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
                
                [maxa, ind] = max( Accuracy*100 )
                maxpara = para(ind)
                maxAcc = max( max(Accuracy*100) )
                %                 matname = sprintf(['C:/Users/csjunxu/Desktop/SC/Results/YaleB_' SegmentationMethod '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_mu' num2str(Par.mu) '_lambda' num2str(lambda) '.mat']);
                %                 save(matname,'Accuracy','medallmissrate','missrateTot','avgmissrate','maxAcc');
            end
        end
    end
end




% rmpath( genpath(currentpath) ) ;
