% Face clustering on the Extended Yale B database by using LSR in

% Can-Yi Lu, Hai Min, Zhong-Qiu Zhao, Lin Zhu, De-Shuang Huang and Shuicheng Yan.
% Robust and Efficient Subspace Segmentation via Least Squares Regression,
% European Conference on Computer Vision (ECCV), 2012.

%--------------------------------------------------------------------------
% Copyright @ Can-Yi Lu, 2012
%--------------------------------------------------------------------------


clear ;

load 'C:\Users\csjunxu\Desktop\SC\Datasets\YaleB_Crop.mat'              % load YaleB dataset
writefilepath = 'C:/Users/csjunxu/Desktop/SC/Results/';
% writefilepath = '';

%% Subspace segmentation methods
% SegmentationMethod = 'LSR1' ;     % LSR1 by (16) in our paper
SegmentationMethod = 'LSR2' ;     % LSR2 by (18) in our paper

Repeat = 20; %number of repeations
DR = 1; % perform dimension reduction or not


%% Data YaleB
nSet = [2:1:10];


%% Subspace segmentation
for set = 1:length(nSet)
    n = nSet(set);
    index = Ind{n};
    %% Parameter
    switch n
        case 5
            para = [0.4] * ones(1,20) ;
        case 10
            para = [0.004 ] * ones(1,20) ;
    end
    for lambda = [1e-6 5e-6 1e-5 5e-5 1e-4 5e-4 1e-3 5e-3 1e-2 5e-2 .1 .5 1 5]
        para = lambda * ones(1,Repeat) ;
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
                redDim = n * 6 ;
            end
            %% normalize
            for c = 1 : size(fea,2)
                fea(:,c) = fea(:,c) /norm(fea(:,c)) ;
            end
            
            %% Subspace segmentation
            missrate = zeros(size(index, 1), Repeat) ;
            
            fprintf( 'dimension = %d \n', redDim ) ;
            Yfea = fea(1:redDim,:) ;
            for j = 1 : Repeat
                p = para( j ) ;
                Accuracy(i,j) = SubspaceSegmentation( SegmentationMethod , Yfea , gnd , p ) ;
                missrate(i, j) = 1 - Accuracy(i,j);
                fprintf('%.3f%% \n' , missrate(i, j)*100) ;
            end
            
            missrateTot{n}(i) = mean(missrate(i, :)*100);
            fprintf('Mean Accuracy of %d/%d is %.3f%%.\n ' , i, size(index, 1), missrateTot{n}(i)) ;
        end
        %% output
        avgmissrate(n) = mean(missrateTot{n});
        medmissrate(n) = median(missrateTot{n});
        fprintf('Total mean missrate  is %.3f%%.\n ' , avgmissrate(n)) ;
        matname = sprintf([writefilepath 'YaleB_Crop_' SegmentationMethod '_DR' num2str(redDim) '_lambda' num2str(lambda) '.mat']);
        save(matname,'missrateTot','avgmissrate','medmissrate');
    end
end
