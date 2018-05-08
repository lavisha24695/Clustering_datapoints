function f = fun1(y, theta, w, N)
    for i=1:1:N
       f(i) = theta.*y - y(i)* y.*w(i,:); 
    end    

end