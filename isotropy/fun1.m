function f = fun1(y, theta, w, N)
    %f = zeros(N,1);
    for i=1:1:N
       f(i) = theta(i) - sum(y.*w(i,:)'); 
    end
    for i=1:1:N
       f(i+N) = y(i)*y(i) - 1; 
    end
    

end


%{ 
%for minimisation
function f = fun1(y, theta, w, N)
    %f = zeros(N,1);
    for i=1:1:N
       f = sum(theta.*y - y(i)*y.*w(i,:)'); 
    end    

end
%}