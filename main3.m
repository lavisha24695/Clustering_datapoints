clc;close all;clear all;
addpath isotropy;addpath misc;addpath hierarchy;
%filename = 'data_spiral.txt';
%filename = 'data_crescents.txt';
%filename = 'sample_data.txt';
%filename = 'data_set_1_8_1_gaussian.txt';
filename = 'Aggregation.txt'
fname = sprintf('rich%s.mat',filename);
load(fname)
%data is the normalized rich data, %data_rich is the unnormalized rich data
%data_original is the original data in the original space
no_properties = size(data,2);
lr = 0.000001; %learning rate
w = rand(no_properties,1);
w = softmax(w);
%For each dataset do the following
N = size(data,1);
D = size(data,2);
%neighbors_criterion = 'Fixed_number';
neighbors_criterion = 'Fixed_radius';
% Precompute Neighbors
if strcmp(neighbors_criterion, 'Fixed_number')
    data_neighbors = precompute_nearest_neighbors(data, params.Kmax);
else
    [data_neighbors, nbr_radius] = precompute_neighbors_fixedradius(data_original);        
end

y = data(:,10);%Uniformity isotropy is the 10th dimension
%meany = median(y_orig);%stdy = std(y_orig);y = y_orig > meany - stdy;
alpha = 0.05;
y = y > alpha;
f = figure;
scatter(data_original(:,1), data_original(:,2), 1,y);

T = 100;
weights = zeros(T+1,no_properties);
weights(1,:) = w';
labels = zeros(T, N);
iter = 1;
flow = zeros(iter,1);
loss = zeros(iter,1);
for iter = 1:1:T
    weight=weights(iter,:);
    [yhat, flow(iter), unary, binary, binary2]= energy_minimisation(data, weight,data_neighbors, data_original);
    yhat = logical(yhat);
    labels(iter,:) = yhat';
    f = figure;
    scatter(data_original(:,1), data_original(:,2), 1,yhat);
    loss(iter) = sumsqr(yhat - y);
    fprintf('\nIteration: %d , loss: %d, flow: %d', iter, loss(iter), flow(iter));
    derivative = compute_derivate(weights(iter,:), data, y, yhat, data_neighbors);
    %fprintf('\n Old weights: %d',  weights(iter,:)); fprintf('\n Derivative: %d ', lr*derivative');
    
    w = weights(iter,:) - lr*derivative';
    w2 = (w-min(w))/(max(w)-min(w));
    weights(iter+1,:) = softmax(w2');
    weights(iter+1,:) = softmax(rand(no_properties,1));
end
figure;
plot(loss, 'r')
hold on
plot(flow, 'b')
%{
%Visualising the neighbors in the high dimensional space
for i = 1:1:N
    close all;
    i = 321;
   nbrs = data_neighbors{i};
   plot(data_original(:,1), data_original(:,2), 'r.');
   hold on;
   plot(data_original(nbrs,1), data_original(nbrs,2), 'b.');
   hold on;
   plot(data_original(i,1), data_original(i,2), 'g.');
end
%}