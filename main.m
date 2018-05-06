% Script to run pshift clustering.
% Author: Sanketh Shetty (2009-2011)

clc;
close all;
clear all;
addpath isotropy;
addpath misc;
addpath hierarchy;

%filename = 'data_spiral.txt';
filename = 'data_crescents.txt';
%filename = 'sample_data.txt';
%filename = 'data_set_1_8_1_gaussian.txt';
% Example data set:
data = load(filename);

% isotropy_criterion: specify the criterion to use when computing shift
% vectors. Options:
% 'force-sum' : Performs clustering using force-sum criterion
% 'force-sign' : Performs clustering using force-sign criterion
% 'uniformity' : Performs clustering using uniformity criterion
criterion = 'uniformity';

% alpha : specify level of significance of the isotropy test. 
% With the force-sum and force-sign criterion (unpublished), 
% votes are weighted by degree
% of isotropy of a neighborhood. Here a single alpha value determines the
% cutoff point for the voting. If increasing the neighborhood size reduces
% the degree of isotropy below this threshold, voting is stopped.
% e.g. alpha = 0.01;
%lavisha - crescent data, i got perfect labelling with alpha = 0.001
alpha = 0.05;

% win_size: specify a window size for testing.
% Typically set to 25 for force_sign and uniformity.
% For force_sum set it to 10 and the criterion automatically adjusts to the
% number of points used in the test.
win_size = 25;

% algorithm: specify the algorithm to use with the shift vectors to compute
% the final clustering. Options:
% connected-components
% density-hierarchy
% isotropy-hierarchy
algorithm = 'connected-components';

options.min_cluster_size = 20;
options.merge_outliers = 1;

[pvals, data_neighbors, data_dimensions] ...
    = isotropy_clustering(data, criterion, win_size);

if strcmp(algorithm, 'isotropy-hierarchy')
    labels_hierarchy = isotropy_hierarchy(data_neighbors, pvals, win_size);
else
    isborder = pvals < alpha;
    labels = cluster_labeling_mst(isborder,...
            data_neighbors, data, options.min_cluster_size);
end

if strcmp(algorithm, 'density-hierarchy')
    if length(unique(labels)) == 1
        fprintf('Only one cluster detected\n');
        labels_hierarchy{1} = labels;
    else
        [labels_hierarchy, persistence] = ...
        build_hierarchy(labels, data_neighbors, data_dimensions, win_size);  
    end
    
end


% Enjoy!
% Questions & Feedback: sanketh.shetty@gmail.com
figure;
plot3(data(:,1), data(:,2), pvals(:), '.')
figure;
plot(data(isborder,1), data(isborder,2),'.', 'color','r')
hold on
plot(data(~isborder,1), data(~isborder,2),'.', 'color','b')
%
figure;
uniquelabels = unique(labels);
for i=1:size(uniquelabels,1)
    indices = find(~(labels-uniquelabels(i)));
    clusterpts = data(indices,:);
    %figure;
    hold on
    plot(clusterpts(:,1), clusterpts(:,2), '.','color',rand(1,3));
    %axis([-50 20 -15 35])
end
%
%{
for i = 1:size(persistence)
    newlabels = labels_hierarchy{i};
    i
    uniquelabels = unique(newlabels);
    size(uniquelabels,1)
    figure;
    for i=1:size(uniquelabels,1)
        indices = find(~(labels-uniquelabels(i)));
        clusterpts = data(indices,:);
        %figure;
        hold on
        plot(clusterpts(:,1), clusterpts(:,2), '.','color',rand(1,3));
        %axis([-50 20 -15 35])
    end
end
%}