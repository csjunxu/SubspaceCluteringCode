clear;

% %% YaleB = 'YaleB.mat';
% load 'YaleBa.mat';

%% USPS
load 'USPS.mat';
fea = fea';

% %% MNIST
% load 'MNIST.mat';
% fea = fea';



nSet = max(gnd);

Y = cell(nSet, 1);
S = cell(nSet, 1);
for i = 1:nSet
    Y{i, 1} = fea(:, gnd==i);
    S{i, 1} = gnd(find(gnd==i));
end

save USPSa.mat Y S USPSInd;
    