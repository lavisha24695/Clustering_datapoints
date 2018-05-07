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
fname = sprintf('rich%s.mat',filename);
load(fname)
%data is the normalized rich data
%data_rich is the unnormalized rich data
%data_original is the original data in the original space
no_properties = size(data,2);
T = 100;
alpha = 0.1; %learning rate
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
    [data_neighbors, nbr_radius] = precompute_neighbors_fixedradius(data);        
end

%Uniformity isotropy is the 10th dimension
y_orig = data(:,10);
meany = mean(y_orig);
y = y_orig > meany;
weights = zeros(T+1,no_properties);
weights(1,:) = w';
for iter = 1:1:T
    %Minimise energy to find yhat (n X 1)
    %yhat = energy_minimisation(data, weights(iter,:),data_neighbors);
    yhat = data(:,9);
    yhat = yhat > mean(yhat);
    loss = sumsqr(yhat - y)/N;
    fprintf('Iteration: %d , loss: %d', iter, loss);
    %2. Compute weights w using grad-descent
    derivative = compute_derivate(weights(iter,:), data, y, yhat, data_neighbors);
    w = weights(iter,:) - alpha*derivative;
    weights(iter+1,:) = softmax(w);
end