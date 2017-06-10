clear ;
% close all;
currentpath = cd ;
AddedPath = genpath( currentpath ) ;
addpath( AddedPath ) ;
fprintf('\n\n**************************************   %s   *************************************\n' , datestr(now) );
fprintf( [ mfilename(currentpath) ' Begins.\n' ] ) ;
fprintf( [ mfilename(currentpath) ' is going, please wait...\n' ] ) ;

%% reduced dimension
ProjRank = 12 ;
datadir = 'C:/Users/csjunxu/Desktop/SC/Datasets/Hopkins155/';
seqs = dir(datadir);
% Get rid of the two directories: "." and ".."
seq3 = seqs(3:end);
% Save the data loaded in struct "data"
data = struct('ProjX', {}, 'name',{}, 'ids',{});


for i=1:length(seq3)
    fname = seq3(i).name;
    fdir = [datadir '/' fname];
    if isdir(fdir)
        datai = load([fdir '/' fname '_truth.mat']);
        id = length(data)+1;
        % the true group numbers
        data(id).ids = datai.s;
        % file name
        data(id).name = lower(fname);
        % X is the motion sequence
        X = reshape(permute(datai.x(1:2,:,:),[1 3 2]), 2*datai.frames, datai.points);
        
        % PCA projection
        [ eigvector , eigvalue ] = PCA( X ) ;
        ProjX = eigvector(:,1:ProjRank)' * X ;
        data(id).ProjX = [ProjX ; ones(1,size(ProjX,2)) ] ;
    end
end
clear seq3;


%% Subspace segmentation methods
% SegmentationMethod = 'LSR1' ;     % LSR1 by (16) in our paper
% SegmentationMethod = 'LSR2' ;     % LSR2 by (18) in our paper

SegmentationMethod = 'LSRd0po' ;
% SegmentationMethod = 'LSRpo' ;

%% Parameter
switch SegmentationMethod
    case 'LSR1'
        lambda = 4.8*1e-3 ;
    case 'LSR2'
        lambda = 4.6*1e-3 ;
    case 'LSRd0po'
        lambda = 4.8*1e-3 ;
    case 'LSRpo'
        lambda = 4.6*1e-3 ;
end


for mu = [1]
    Par.mu = mu;
    for maxIter = [200 500 1000]
        Par.maxIter = maxIter;
        for rho = [0.01:0.01:0.1]
            Par.rho = rho;
            for lambda = [1e-5:1e-5:5e-5]
                Par.lambda = lambda;
                maxNumGroup = 5;
                for i = 1:maxNumGroup
                    num(i) = 0;
                end
                %%
                errs = zeros(length(data),1);
                for i = 1 : length(data)
                    ProjX = data(i).ProjX ;
                    gnd = data(i).ids' ;
                    K = length( unique( gnd ) ) ;
                    n = max(gnd);
                    switch SegmentationMethod
                        case 'LSRd0po'
                            C = LSRd0po( ProjX , Par ) ;
                        case 'LSRpo'
                            C = LSRpo( ProjX , Par ) ;
                    end
                    nCluster = length( unique( gnd ) ) ;
                    Z = ( abs(C) + abs(C') ) / 2 ;
                    idx = clu_ncut(Z,nCluster) ;
                    accuracy = compacc(idx,gnd) ;
                    missrate = 1-accuracy;
                    num(n) = num(n) + 1;
                    missrateTot{n}(num(n)) = missrate;
                    fprintf('seq %d\t %f\n', i , missrate ) ;
                end
                fprintf('\n') ;
                
                L = [2 3];
                allmissrate = [];
                for i = 1:length(L)
                    j = L(i);
                    avgmissrate(j) = mean(missrateTot{j});
                    medmissrate(j) = median(missrateTot{j});
                    allmissrate = [allmissrate missrateTot{j}];
                end
                avgallmissrate = sum(allmissrate)/length(allmissrate);
                medallmissrate = median(allmissrate);
                matname = sprintf(['C:/Users/csjunxu/Desktop/SC/Results/' SegmentationMethod '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_mu' num2str(Par.mu) '_lambda' num2str(lambda) '.mat']);
                save(matname,'avgallmissrate','medallmissrate','missrateTot','avgmissrate','medmissrate');
            end
        end
    end
end




