function derivative = compute_derivate(weight, data, y, yhat, data_neighbors)
    D = size(weight,2);
    N = size(data,1);
    derivative = zeros(D,1);
    
    %Compute wij
    w = zeros(N,N);
    for i = 1:1:N
        nbrs = data_neighbors{i};
        w(i,i) = 1;
        for j = 1:1:size(nbrs,2)
            w(i,nbrs(j)) = exp(-1* sum(weight.*(data(i,:) - data(nbrs(j),:)).*(data(i,:) - data(nbrs(j),:))));
        end
    end
    fprintf('\nComputed w\n');
   %Computing thetai
   p = zeros(N,1);
   theta = zeros(N,1);
   for i = 1:1:N
       p(i) = sum(data(i,:).*weight);
   end
   [pmax, maxi] = max(p);
   [pmin, mini] = min(p);
   a = maxi;
   b = mini;
  % fprintf('Pmax: %d, Pmin: %d\n',pmax, pmin);
   for i = 1:1:N
       theta(i) = 0.5*(pmax + pmin - 2*p(i))/(pmax - pmin);
   end
   dldyhat = 2*(yhat-y);
   
   %Compute dEdyhat
   dEdyhat = zeros(N,1);
   for i = 1:1:N
       nbrs =data_neighbors{i};
       dEdyhat(i) = theta(i) - sum(yhat(nbrs)'.*w(i, nbrs));
   end
   
   fprintf('Compute dEdweightd\n');
   dEdweightd = zeros(D,1);
   for d = 1:1:D
       sum1 = 0;
       for i = 1:1:N
           nbrs = data_neighbors{i};
           part2 = sum(2*yhat(i).*yhat(nbrs)'.*w(i,nbrs).*(repmat(data(i,d), size(nbrs)) - data(nbrs,d)'));
           dthetaidweightd = pmax*data(b,d) - pmin*data(a,d) + p(i)*(data(a,d) - data(b,d)) - data(i,d)*(pmax-pmin);
           dthetaidweightd = dthetaidweightd/((pmax-pmin)*(pmax-pmin));
           part1 = yhat(i)*dthetaidweightd; 
           sum1 = sum1 + part1 + part2;
       end
       dEdweightd(d) = sum1;
   end
   dEdyhat(find(dEdyhat==0)) = 1;
   for d = 1:1:D
       derivative(d) = dEdweightd(d)*sum(dldyhat./dEdyhat);
   end
end