function gradient = compute_fungradient(i, data, data_nbrs, fun)

d =  size(data,2);
n = size(data_nbrs,2)+1;
ptwithnbrs = zeros(n,d+1);
ptwithnbrs(1:n-1,1:d) = data(data_nbrs,:);
ptwithnbrs(1:n-1,end) = fun(data_nbrs);
ptwithnbrs(end,1:d) = data(i);
ptwithnbrs(end,end) = fun(i);

grad_dims = zeros(d,1);
for j = 1:1:d
    A1 = sortrows(ptwithnbrs,j);
    A2 = zeros(size(A1));
    A2(1:n-1,:) = A1(2:n,:);
    diff = A2-A1;
    grads = diff(:,d+1)./diff(:,j);
    grads(isnan(grads)) = 0;
    grads(isinf(grads)) = 0;
    
    grad_dims(j) = mean(grads);
    %grad_dims(j) = nansum(grads)/size(grads,1);
end
gradient = mean(grad_dims);
if isnan(gradient)
    fprintf('%d, nan encountered\n',i);
end
end