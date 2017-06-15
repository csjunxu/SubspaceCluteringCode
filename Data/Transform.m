clear;

load 'YaleB.mat';

nSet = max(gnd);

Y = cell(nSet, 1);
S = cell(nSet, 1);
for i = 1:nSet
    Y{i, 1} = fea(:, find(gnd==i));
    S{i, 1} = gnd(find(gnd==i));
end

save YaleB_Crop.mat Y S Ind;
    