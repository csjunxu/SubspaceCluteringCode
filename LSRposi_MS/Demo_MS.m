clear ;

%% reduced dimension
ProjRank = 12 ;
datadir = 'C:/Users/csjunxu/Desktop/SC/Datasets/Hopkins155/';
seqs = dir(datadir);
% Get rid of the two directories: "." and ".."
seq3 = seqs(3:end);
% Save the data loaded in struct "data "
data = struct('ProjX', {}, 'name',{}, 'ids',{});

dataset = 'Hopkins155';

resultdir = 'C:/Users/csjunxu/Desktop/SC/Results/';
 
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
% SegmentationMethod = 'LSR' ;
% SegmentationMethod = 'LSRd0' ;

% SegmentationMethod = 'NNLSR_LSR' ;
% SegmentationMethod = 'NNLSRd0_LSR' ;

SegmentationMethod = 'ANNLSR_LSR' ;
% SegmentationMethod = 'ANNLSRd0_LSR' ;

for mu = [1]
    Par.mu = mu;
    for maxIter = [5 10 15]
        Par.maxIter = maxIter;
        for rho = [.006 .005 0.004]
            Par.rho = rho;
            for lambda = [1 2 3 4 5 10 .1 .2 .3 .4 .5]
                Par.lambda = lambda*10^(-3);
                maxNumGroup = 5;
                for i = 1:maxNumGroup
                    num(i) = 0;s
                end
                %%
                errs = zeros(length(data),1);
                for i = 1 : length(data)
                    ProjX = data(i).ProjX ;
                    gnd = data(i).ids' ;
                    K = length( unique( gnd ) ) ;
                    n = max(gnd);
                    switch SegmentationMethod
                        case 'LSR'
                            C = LSR( ProjX , Par ) ;
                        case 'LSRd0'
                            C = LSRd0( ProjX , Par ) ; % solved by ADMM
                            % C = LSR1( ProjX , Par.lambda ) ; % proposed by Lu
                        case 'NNLSR_LSR'                   % non-negative
                            C = NNLSR( ProjX , Par ) ;
                        case 'NNLSRd0_LSR'               % non-negative, sum to 1
                            C = NNLSRd0( ProjX , Par ) ;
                        case 'ANNLSR_LSR'                 % affine, non-negative
                            C = ANNLSR( ProjX , Par ) ;
                        case 'ANNLSRd0_LSR'             % affine, non-negative, sum to 1
                            C = ANNLSRd0( ProjX , Par ) ;
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
                matname = sprintf([resultdir dataset '_' SegmentationMethod '_maxIter' num2str(Par.maxIter) '_rho' num2str(Par.rho) '_mu' num2str(Par.mu) '_lambda' num2str(Par.lambda) '.mat']);
                save(matname,'avgallmissrate','medallmissrate','missrateTot','avgmissrate','medmissrate');
            end
        end
    end
end




