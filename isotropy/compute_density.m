function density = compute_density(i, data, data_nbrs, nbr_radius)

%nbrs = data(data_nbrs,:);
n = size(data_nbrs,2);
d = size(data,2);
volume = 1;
constant = 1;
for i = 1:1:d
    volume = volume*nbr_radius;
end
density = constant*n/volume;

%The code below computes the median distance of the given point i from all
%its neighbors
%{
c = 1;
distance = zeros(n,1);
for j = 1:1:n
    vect = data(i,:) - nbrs(j,:);
    distance(j) = sqrt(sumsqr(vect));
end
dc = median(distance);
if dc==0
    density = 0;
else
    density = c/(dc*dc);
end    
%}

end