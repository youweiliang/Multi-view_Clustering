function [distance_matrix] = make_distance_matrix_single(data, metric)
% Construct the weight matrix (a graph) for data.
% Inputs:
%   data - a data feature matrix with each row being an instance (data point), each column representing a feature
%   metric - the metric used to measure the distance between two instances
% Outputs:
%   affinity_matrix - a matrix representing the similarity between data points
%   distance_matrix - a matrix representing the distance (dissimilarity) between data points

if strcmp(metric, 'cosine')
    distance_matrix = pdist2_fast(data, data, 'cosine');
    distance_matrix = distance_matrix.^2;
elseif strcmp(metric, 'original')
    distance_matrix = data;
elseif strcmp(metric, 'squaredeuclidean')
    distance_matrix = pdist2(data, data, 'squaredeuclidean');
else
    error('unknown metric')
end 

end