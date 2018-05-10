clc; close all; clear all;
addpath isotropy; addpath misc; addpath hierarchy;

filename = 'Aggregation.txt'
%filename = 'data_spiral.txt';
%filename = 'data_crescents.txt';
%filename = 'sample_data.txt';
%filename = 'data_set_1_8_1_gaussian.txt';
data = load(filename);
data = data(:,1:2);

% 'force-sum' : Performs clustering using force-sum criterion
% 'force-sign' : Performs clustering using force-sign criterion
% 'uniformity' : Performs clustering using uniformity criterion
criterion = 'force-sign';
%lavisha - crescent data, i got perfect labelling with alpha = 0.001
alpha = 0.05;
win_size = 25;
algorithm = 'connected-components';

options.min_cluster_size = 20;
options.merge_outliers = 1;

[pvals, data_neighbors, data_dimensions] = isotropy_clustering(data, criterion, win_size);

%fname = sprintf('rich%s.mat',filename);
%load(fname)
%pvals = (data(:,4));
%alpha = mean(pvals);
%data = data_original;

alpha = 0.05;
isborder = pvals > alpha;
figure;
plot(data(isborder,1), data(isborder,2),'.', 'color','r')
hold on
plot(data(~isborder,1), data(~isborder,2),'.', 'color','b')

labels = cluster_labeling_mst(isborder,data_neighbors, data, options.min_cluster_size);

figure;
plot3(data(:,1), data(:,2), pvals(:), '.')

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
