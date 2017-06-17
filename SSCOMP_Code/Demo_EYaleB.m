% This code tests the performance of SSC-OMP on the EYaleB face database.
% The code generates results in Table 2 of the paper
% C. You, D. Robinson, R. Vidal, Scalable Sparse Subspace Clustering by
% Orthogonal Matching Pursuit, CVPR 2016.

% In this code, we apply SSC-OMP to the face images of a randomly picked
% n \in \{2, 10, 20, 30, 38} subjects in the Extended Yale B database. Each
% subject has 64 images under different illumination conditions.

% Instructions for running the code:
% - Download code for computing clustering accuracy. Go to
% http://www.cad.zju.edu.cn/home/dengcai/Data/Clustering.html and download
% bestMap.m and Hungarian.m. Alternatively, you can use your own function
% by redefining the function evalAccuracy.m
% - Run. You can modify the parameter "nCluster" below to run for different
% number of subjects.

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com

%% Settings
% setup
nCluster = 2; % number of subjects.

%% Load data
% Data is preprocessed and saved in the .mat file.
% EYALEB_DATA is a D by N matrix. Each column is a face image and N =
% 38 subjects * 64 images/subject = 2414. Each image is downsampled from
% 192*168 to D = 48*42 = 2016.
% EYALEB_LABEL is a 1 by N vector. Each entry is the label for the
% corresponding column in EYALEB_DATA.

load 'C:\Users\csjunxu\Desktop\SC\Datasets\YaleB_Crop.mat'              % load YaleB dataset
dim = 6;


% dimension reduction
reduceDimension = @(data) dimReduction_PCA(data, dim*nCluster); % 0
% normalization
normalizeColumn = @(data) cnormalize_inplace(data);
% representation
buildRepresentation = @(data) OMP_mat_func(data, 5, 1e-6); % second parameter is sparsity
% spectral clustering
genLabel = @(affinity, nCluster) SpectralClustering(affinity, nCluster, 'Eig_Solver', 'eigs');


%% Clustering
Repeat = 1; %number of repeations
results = zeros(nExperiment, 6); %results

for nSet = [2 3 5 8 10];
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
        for j = 1 : Repeat         
            N = length(gnd);
            % clustering
            tic;
            %     fprintf('Dimension reduction...\n')
            fea = reduceDimension(fea);
            % normalization
            %     fprintf('Normalization...\n')
            fea = normalizeColumn(fea);
            % generate representation
            %     fprintf('Representation...\n')
            R = buildRepresentation(fea);
            % generate affinity
            %     fprintf('Affinity...\n')
            R(1:N+1:end) = 0;
            % R = cnormalize(R, Inf);
            A = abs(R) + abs(R)';
            % generate label
            %     fprintf('Generate label...\n')
            groups = genLabel(A, nCluster);
            time = toc;
            
            % Evaluation
            perc = evalSSR_perc( R, gnd );
            ssr = evalSSR_error( R, gnd );
            conn = evalConn( A, gnd);
            accr  = evalAccuracy(gnd, groups);
            % output
            dataformat = '%d-th experiment: perc = %f, ssr = %f, conn = %f, accr = %f, time = %f\n';
            dataValue = [iExperiment, perc, ssr, conn, accr, time];
            fprintf(dataformat, dataValue);
            % record
            results(iExperiment, :) = dataValue;
        end
        % output
        dataValue = mean(results, 1);
        fprintf('\nAverage: perc = %f, ssr = %f, conn = %f, accr = %f, time = %f\n', dataValue(2:end));
    end
end
