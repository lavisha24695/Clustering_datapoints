function shiftvector = compute_shiftvector(i, data, data_nbrs)

n = size(data_nbrs,2);
d = size(data,2);
sum1 = zeros(1,d);

for j=1:1:n
   sum1 = sum1 + data(i,:) - data(data_nbrs(j),:); 
end

sum1 = sum1./n;
shiftvector = sqrt(sumsqr(sum1));

end