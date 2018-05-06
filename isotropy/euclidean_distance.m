%Function to compute Euclidean distance between point P and data C
function D = euclidean_distance(P,C)


X1 = P(:,ones(size(C,2),1));

B = (C-X1).*(C-X1);
D = sum(B,1);

D = sqrt(D);

end
