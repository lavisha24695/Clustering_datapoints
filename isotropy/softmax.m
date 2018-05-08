function weightsnew = softmax(w2)
    D = size(w2,1);
    sum1 = 0;
    for i = 1:1:D
        sum1 = sum1 + exp(w2(i));
    end
    weightsnew = exp(w2)./sum1;
end