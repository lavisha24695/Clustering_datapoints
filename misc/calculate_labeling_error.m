%% Function to compute the clustering error.
%  Author: Sanketh Shetty (September 2010)

function [unlabeled,average_per_cluster_impurity,total_cluster_impurity]...
      = calculate_labeling_error(expected_labels,data_labels)

    unlabeled = sum(data_labels == 0);
   
    labels = unique(data_labels);
    
    if labels(1) == 0
        labels(1) = [ ];
    end
    
    average_per_cluster_impurity = 0;
    total_cluster_impurity = 0;
    
    for i=1:length(labels)
        
        ids = find(data_labels==labels(i));
        l2 = unique(expected_labels(ids));
        
        % No error.
        if length(l2) == 1
            continue;
        end
        
        N = hist(expected_labels(ids),l2);
        [max_set,max_l]=max(N);
        
        average_per_cluster_impurity = average_per_cluster_impurity + (sum(N)-max_set)/sum(N);
        total_cluster_impurity = total_cluster_impurity + (sum(N)-max_set);
    end
    
    average_per_cluster_impurity = average_per_cluster_impurity/length(labels);
    total_cluster_impurity = total_cluster_impurity/length(expected_labels);
end
