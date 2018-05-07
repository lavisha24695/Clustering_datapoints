function derivative = compute_derivate(weights, data, y, yhat, data_neighbors)
    D = size(weights,1);
    N = size(data,1);
    derivative = zeros(D,1);
    %Compute wij
    w = zeros(N,N);
    for i = 1:1:N
        for j = 1:1:N
            w(i,j) = exp(-1* sum(weights.*(data(i,:) - data(j,:)).*(data(i,:) - data(j,:))));
            w(j,i) = exp(-1* sum(weights.*(data(i,:) - data(j,:)).*(data(i,:) - data(j,:))));

        end
    end
   %Computing thetai
   p = zeros(N,1);
   theta = zeros(N,1);
   for i = 1:1:N
       p(i) = sum(data(i,:).*weights);
   end
   [pmax, maxi] = max(p);
   [pmin, mini] = min(p);
   a = maxi;
   b = mini;
   fprintf(maxv, minv);
   for i = 1:1:N
       theta(i) = 0.5*(p(maxi) + p(mini) - 2*p(i))/(p(maxi) - p(mini));
   end
   dldyhat = 2*(yhat-y);
   %Compute dEdyhat
   dEdyhat = zeros(N,1);
   for i = 1:1:N
       nbrs =data_neighbors{i};
       dEdyhat(i) = theta(i) - sum(yhat(nbrs).*w(i, nbrs));
   end
   %Compute dwijdweightd
  
   %Compute dEdweightd
   dEdweightd = zeros(D,1);
   for d = 1:1:D
       sum1 = 0;
       for i = 1:1:N
           nbrs =data_neighbors{i};
           part2 = sum(2*yhat(i)*yhat(nbrs).*w(i,nbrs).*( repmat(data(i,d), size(nbrs),1) - data(nbrs,d)));
           dthetaidweightd = 2*pmax*data(b,d) - 2*pmin*data(a,d) + 2*p(i)*(data(a,d) - data(b,d)) - 2*data(i,d)*(pmax-pmin);
           dthetaidweightd = dthetaidweightd/(2*(pmax-pmin)*(pmax-pmin));
           part1 = yhat(i)*dthetaidweightd; 
           sum1 = sum1 + part1 + part2;
       end
       dEdweightd(d) = sum1;
   end
   
   for d = 1:1:D
       derivative(d) = dEdweightd(d)*sum(dldyhat./dEdyhat);
   end
end