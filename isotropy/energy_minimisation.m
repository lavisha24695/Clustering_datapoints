function [labels,flow] = energy_minimisation(data, weight, data_neighbors)
   
    D = size(weight,2);
    N = size(data,1);
    % y = zeros(n,1);
    w = zeros(N,N);
    maxnbr = 0;
    for i = 1:1:N
        nbrs = data_neighbors{i};
        %w(i,i) = 1;
        for j = 1:1:size(nbrs,2)
            w(i,nbrs(j)) = exp(-1* sum(weight.*(data(i,:) - data(nbrs(j),:)).*(data(i,:) - data(nbrs(j),:))));
        end
        if(size(nbrs,2)>maxnbr)
            maxnbr = size(nbrs,2);
        end
    end
    p = zeros(N,1);
    theta = zeros(N,1);
    for i = 1:1:N
       p(i) = sum(data(i,:).*weight);
    end
    [pmax, maxi] = max(p);
    [pmin, mini] = min(p);
    theta1 = zeros(N,2);
    %theta2 = zeros(N,1);
    fprintf('Pmax: %d, Pmin: %d\n',pmax, pmin);
    for i = 1:1:N
       theta1(i,1) = 0.5*(pmax - p(i))/(pmax-pmin);
       theta1(i,2) = 0.5*(p(i) - pmin)/(pmax-pmin);
       theta(i) = 0.5*(pmax + pmin - 2*p(i))/(pmax - pmin);
    end
%    w2 = w;
 %   theta1 = 10*theta1;
    A = sparse(w);
    T = sparse(theta1);
    [flow, labels] = maxflow(A,T);
    %{
        E = edges4connected(height,width);
        V = abs(m(E(:,1))-m(E(:,2)))+eps;
        A = sparse(E(:,1),E(:,2),V,N,N,4*N);
        % terminal weights
        % connect source to leftmost column.
        % connect rightmost column to target.
        T = sparse([1:height;N-height+1:N]',[ones(height,1);ones(height,1)*2],ones(2*height,1)*9e9);
        disp('calculating maximum flow');
        %A is N*N
        %T is N*2
        [flow,labels] = maxflow(A,T);
    %}
    %y0 = randi([0 1], N,1);
    %f = @(y)fun1(y, theta, w, N);
    %[y, fval] = fminunc(f,y0, @mycon);
    %[y, fval] = fsolve(f,y0);
end