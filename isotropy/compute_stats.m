function [avg, stddev] = compute_stats(i, data, data_nbrs)

n = size(data_nbrs,2)+1;
ptwithnbrs(1:n-1,:) = data(data_nbrs,:);
ptwithnbrs(n,:) = data(i);

[ndata, dimx] = size(ptwithnbrs);
[ncentres, dimc] = size(ptwithnbrs);
n2 = (ones(ncentres, 1) * sum((ptwithnbrs.^2)', 1))' + ...
  ones(ndata, 1) * sum((ptwithnbrs.^2)',1) - ...
  2.*(ptwithnbrs*(ptwithnbrs'));
if any(any(n2<0))
  n2(n2<0) = 0;
end
dist = sqrt(n2);

%n1 = dist/(ndata-1);
%nbr_radius = sum(sum(n1))/ndata;
avg = mean(dist(:));
stddev = std(dist(:));
%density = 1/nbr_radius;

end