% Function to output status of loops


function diagnostics(i)

    if mod(i,100) == 0
        fprintf('.');
    end

    if mod(i,1000) == 0
        fprintf('%d\n',i);
    end
    
end
