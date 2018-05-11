clc; close all; clear all;
addpath isotropy; addpath misc; addpath hierarchy;

%filename = 'Aggregation.txt'
%filename = 'data_spiral.txt';
filename = 'data_crescents.txt';
%filename = 'sample_data.txt';
%filename = 'data_set_1_8_1_gaussian.txt';
data = load(filename);
data = data(:,1:2);

params.Kmax = min(200, size(data,1)); % maximum window size for testing
data_neighbors = precompute_nearest_neighbors(data, params.Kmax);

fname = sprintf('labels%s.mat',filename);
load(fname)
%{
% 'force-sum' : Performs clustering using force-sum criterion
% 'force-sign' : Performs clustering using force-sign criterion
% 'uniformity' : Performs clustering using uniformity criterion
criterion = 'force-sign';
win_size = 25;
[pvals, data_neighbors, data_dimensions] = isotropy_clustering(data, criterion, win_size);
%fname = sprintf('rich%s.mat',filename);
%load(fname)
%pvals = (data(:,4));
figure;
plot3(data(:,1), data(:,2), pvals(:), '.')
%alpha = mean(pvals);
%data = data_original;
alpha = 0.05;
isborder = pvals > alpha;
%}
f = figure;
plot(data(isborder,1), data(isborder,2),'.', 'color','r')
hold on
plot(data(~isborder,1), data(~isborder,2),'.', 'color','b')
saveas(f, 'boundary.png') 
options.min_cluster_size = 20;
labels = cluster_labeling_mst(isborder,data_neighbors, data, options.min_cluster_size);
f1 = figure;
uniquelabels = unique(labels);
for i=1:size(uniquelabels,1)
    indices = find(~(labels-uniquelabels(i)));
    clusterpts = data(indices,:);
    hold on
    plot(clusterpts(:,1), clusterpts(:,2), '.','color',rand(1,3));
end
saveas(f1, 'cluster.png')