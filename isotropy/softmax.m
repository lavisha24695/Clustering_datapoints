function weightsnew = softmax(weights)
    D = size(weights,1);
    sum1 = 0;
    for i = 1:1:D
        sum1 = sum1 + exp(weights(i));
    end
    weightsnew = weights./sum1;
end