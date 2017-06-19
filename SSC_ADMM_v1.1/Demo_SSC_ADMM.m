
clear;

load 'C:\Users\csjunxu\Desktop\SC\Datasets\YaleB_Crop.mat';
% load 'C:\Users\csjunxu\Desktop\SC\Datasets\USPS_Crop.mat'   % load USPS dataset

writefilepath = 'C:/Users/csjunxu/Desktop/SC/Results/';

%% Subspace segmentation methods
SegmentationMethod = 'SSC_ADMM' ;
DR = 1; % dimension reduction
dim = 6;
%% Subspace segmentation
for alpha = [4 3 2 1 0.8:-0.2:0.2 0.1 0.05 0.01]
    for nSet = [2 3 5 8 10]
        n = nSet;
        index = Ind{n};
        for i = 1:size(index,1)
            fea = [];
            gnd = [];
            for p = 1:n
                fea = [fea Y{index(i, p), 1}];
                gnd= [gnd p*ones(1, length(S{index(i, p)}))];
            end
            [D,N] = size(fea);
            
            redDim = size(fea, 1);
            if DR == 1
                %% PCA Projection
                [ eigvector , eigvalue ] = PCA( fea ) ;
                maxDim = length(eigvalue);
                fea = eigvector' * fea ;
                redDim = nSet * dim ;
            end
            % normalize
            for c = 1 : size(fea,2)
                fea(:,c) = fea(:,c) /norm(fea(:,c)) ;
            end
            Yfea = fea(1:redDim, :) ;
            r = 0; affine = false; outlier = true; rho = 1;
            [missrate,C] = SSC(Yfea,r,affine,alpha,outlier,rho,gnd);
            
            fprintf('missrate: %2.4f\n', missrate);
            missrateTot{n}(i) = missrate;
        end
        avgmissrate(n) = mean(missrateTot{n});
        medmissrate(n) = median(missrateTot{n});
        matname = sprintf([writefilepath 'YaleB_' SegmentationMethod '_DR' num2str(DR) '_dim' num2str(dim) '_alpha' num2str(alpha) '.mat']);
        save(matname,'missrateTot','avgmissrate','medmissrate');
    end
    matname = sprintf([writefilepath 'YaleB_' SegmentationMethod '_DR' num2str(DR) '_dim' num2str(dim) '_alpha' num2str(alpha) '.mat']);
    save(matname,'missrateTot','avgmissrate','medmissrate');
end



