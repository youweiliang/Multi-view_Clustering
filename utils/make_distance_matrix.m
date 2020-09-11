function WW = make_distance_matrix(fea, metric)
if nargin < 2
    metric = 'squaredeuclidean';
end
v = length(fea);
WW = cell(v, 1);

for i=1:v
    WW{i} = make_distance_matrix_single(fea{i}, metric);
end