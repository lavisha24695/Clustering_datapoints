function labels_hierarchy = isotropy_hierarchy(data_neighbors, pvals, win_size)

    num_pts = length(pvals);
    isactive = zeros(num_pts,1);
    [~,order] = sort(pvals,'descend');

    labels = zeros(num_pts,1);
    
    labels_hierarchy{1} = labels; 
    next_level = 2;
    
    for i = 1:num_pts
        next_point = order(i);
        isactive(next_point) = 1;
        labels(next_point) = next_level;
        
        for j = 1:win_size
            curr_nbr = data_neighbors(next_point).order(j);
            if isactive(curr_nbr) == 1 && sum(data_neighbors(curr_nbr).order == next_point) > 0 &&labels(curr_nbr) ~= next_level
                labels(labels==labels(curr_nbr)) = next_level;
            end
        end
        labels_hierarchy{next_level} = labels;
        next_level = next_level + 1;
    end
end