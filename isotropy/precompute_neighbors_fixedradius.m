function [data_neighbors, nbr_radius] = precompute_neighbors_fixedradius(data)
%Author: Lavisha Aggarwal
%This function finds neighbors of all points within a fixed radius around
%the point. The fixed radius is proportional to the mean distance between
%all points in the dataset
fprintf('Computing fixed radius neighborhoods around point\n');
[ndata, dimx] = size(data);
[ncentres, dimc] = size(data);
%c = 0.1 worked best for main2
%c = 0.5 might be working best for main3
c = 0.1;
n2 = (ones(ncentres, 1) * sum((data.^2)', 1))' + ...
  ones(ndata, 1) * sum((data.^2)',1) - ...
  2.*(data*(data'));
% Rounding errors occasionally cause negative entries in n2
if any(any(n2<0))
  n2(n2<0) = 0;
end
dist = sqrt(n2);
n1 = dist/(ndata-1);
nbr_radius = c*sum(sum(n1))/ndata;
fprintf('Fixed Radius: %d \n', nbr_radius);
data_neighbors = cell(ndata,1);
for i=1:1:ndata
   data_neighbors{i} = find(dist(i,:)<nbr_radius & dist(i,:));
end
fprintf('Done\n');

%{
%Visualising the neighborhood around a point
pt = 828;
plot(data(data_neighbors{pt},1), data(data_neighbors{pt},2), 'r.');
hold on
plot(data(pt,1), data(pt,2), 'k.')
hold on
indices = linspace(1,size(data,1), size(data,1));
nonnbr_index = setdiff(indices, data_neighbors\{pt});
plot(data(nonnbr_index,1), data(nonnbr_index,2), 'b.')
%}
end