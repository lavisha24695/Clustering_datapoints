% Function to filter out small clusters and merge them with larger ones

function label = filter_smaller_clusters(data_labels,data_neighbors,min_cluster_size,merge_outliers)


%Seach the clusters for small clusters and remove them. This is
%optional.

start_size = length(unique(data_labels));
total_points = 0;

    label = data_labels;
    
    fprintf('Removing small clusters.\n');
    uniq_labels = unique(label);

    for l=1:length(uniq_labels)
        
        
        if length(find(label == uniq_labels(l))) < min_cluster_size    
            total_points = total_points+length(find(label==uniq_labels(l)));
            label(find(label == uniq_labels(l))) = 0;
        end
        
    end
    
    find_label = find(label == 0);
    
    % If instructions are to not merge outliers the return
    if merge_outliers==0
        return;
    end
    
    change = -1;

    while change~=0
      
      s1 = length(find_label);
      
      label2 = label;
      
      for i=1:length(find_label)

          j=find_label(i);
          %interior_pts1 = find(isborder(data_neighbors(j).order)==0);
          interior_pts = find(label(data_neighbors{j,1})~=0);
                   
          if isempty(interior_pts)
              continue;
          end
          
          label_list =label(data_neighbors{j,1});

          %Find cluster with maximum vote
          unique_labels = unique(label_list);

          %Outlier
          if length(unique_labels) == 0
            continue;
          end
          
          if unique_labels(1) == 0
              unique_labels = unique_labels(2:end);
          end

          if length(unique_labels) == 0
            continue;
          end

          if length(unique_labels) > 1
            N = hist(label_list,unique_labels);
            [max_set,max_l]=max(N);
            label2(j)=unique_labels(max_l);

          else
              label2(j)=unique_labels(1);
          end
      end
      
      label = label2;
      
      find_label = find(label==0);
      s2 = length(find_label);
      change = s2-s1;

    end
    
    end_size = length(unique(label));
    
    fprintf('Removed %d Clusters containing %d points\n',start_size-end_size,total_points);
    
 end