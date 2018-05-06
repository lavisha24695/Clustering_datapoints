%Function for cluster labeling

%November 2008: Modified so that data_neighbors contains exactly the amount
%of information needed
function data_labels = cluster_labeling2(isborder,data_neighbors,filter_small,merge_noise)

    label = 1:length(isborder);
    border_stop = zeros(size(isborder));
    
    next_l = 1;
    interior = find(isborder == 0);
    
%% Connected Components on Cluster Interiors

fprintf('\nFinding Core Neighbors\n');
    for int_pt=1:length(interior)

        i=interior(int_pt);

        %In the consistency neighborhood of the point look for data_neighbors that
        %are mutual neighbors.
        
        if isempty(data_neighbors(i).order)
            border_stop(i) = 1;
            continue;
        end
        
        for j=1:length(data_neighbors(i).order)

            
            %Grow the neighborhood till a border point is reached.
            %These border data_neighbors have been filtered to remove noisy labels.
            if isborder(data_neighbors(i).order(j)) == 1
                break;
            end
            border_stop(i) = j;
        end

    end

    fprintf('... done.\n');
    

    fprintf('\nConnected Components Analysis on Core Neighbors ...\n');
    
    for int_pt=1:length(interior)

        i=interior(int_pt);

        %In the consistency neighborhood of the point look for data_neighbors that
        %are mutual neighbors.
        
        if isempty(data_neighbors(i).order)
            continue;
        end
        
        for j=1:border_stop(i)
            nbr = data_neighbors(i).order(j);
            its_nbrs = data_neighbors(nbr).order(1:border_stop(nbr));
            
            if sum(its_nbrs == i) > 0
                label(label == label(nbr)) = label(i);
            end
        end

    end

    fprintf('... done.\n');
    
    label(isborder) = 0;

    

%% Assign border points to clusters by voting.
% For each border point we find its set of mutual K neighbors. We select
% the subset that is labeled interior. Each member of this set votes for
% ownership of this border point, on behalf of its cluster. The cluster
% with the maximum number of votes claims ownership of the border point. 

    fprintf('Assigning border points to clusters...');
    
    find_label = find(label == 0);
    
    change = -1;

    while change~=0
      
      s1 = length(find_label);
      
      label2 = label;
      
      for i=1:length(find_label)

          j=find_label(i);
          %interior_pts1 = find(isborder(data_neighbors(j).order)==0);
          interior_pts = find(label(data_neighbors(j).order)~=0);
                   
          label_list =[ ];
          
          %Find KNN that are interior and verify that they are mutual
          %neighbors
          
          for k=1:length(interior_pts)
              
              temp_id = data_neighbors(j).order(interior_pts(k));

              %If mutual neighbors cast a vote on behalf of native cluster.
              if ~isempty(find(data_neighbors(temp_id).order==j))
                  %If the point and the border point are mutual neighbors.
                  label_list = [label_list; label(temp_id)];
              end

          end

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
    fprintf('Done\n');
    
    %Merge all zero labels with nearest cluster points
    if merge_noise == 1

        fprintf('Find outliers and merge them with nearest clusters\n');
        find_label = find(label == 0);

        for i=1:length(find_label)
            j=find_label(i);
            
            for k=1:length(data_neighbors(j).order)
            
                if label(data_neighbors(j).order(k))~=0
                label(j)=label(data_neighbors(j).order(k));
                break;
                end
            end
        end
    end

    %% Reassign small clusters as "outliers/border"


%Seach the clusters for small clusters and remove them. This is
%optional.

    if filter_small ~= 0

        fprintf('Removing small clusters.\n');
        uniq_labels = unique(label);

        for l=1:length(uniq_labels)
            if length(find(label == uniq_labels(l))) < filter_small    
                label(find(label == uniq_labels(l))) = 0;
            end
        end
    end

    find_label = find(label == 0);
    
    change = -1;

    while change~=0
      
      s1 = length(find_label);
      
      label2 = label;
      
      for i=1:length(find_label)

          j=find_label(i);
          %interior_pts1 = find(isborder(data_neighbors(j).order)==0);
          interior_pts = find(label(data_neighbors(j).order)~=0);
                   
          label_list =[ ];
          
          %Find KNN that are interior and verify that they are mutual
          %neighbors
          
          for k=1:length(interior_pts)
              
              temp_id = data_neighbors(j).order(interior_pts(k));

              %If mutual neighbors cast a vote on behalf of native cluster.
              if ~isempty(find(data_neighbors(temp_id).order==j))
                  %If the point and the border point are mutual neighbors.
                  label_list = [label_list; label(temp_id)];
              end

          end

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
    fprintf('Done\n');
 
    
    data_labels = label;
end
