function data_labels = cluster_labeling_mst...
    (isborder,data_neighbors,data,filter_small)

    label = (1:length(isborder))';
    border_stop = zeros(size(isborder));
    
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

		if sum(label==0) == length(label)
				data_labels = -ones(1,length(label));
				return;
		end

    %plotClusters([data label'],2,1,0);
    
    label = mst_label(data,label);
    
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
    
		if sum(label==0) > 0
				data_labels = mst_label(data,label);
    else
				data_labels = label;
	  end
end

function label = mst_label(data,label)
    notassigned = find(label==0);
    assigned = find(label~=0);
    
    Q = data(notassigned,:);
    P = data(assigned,:);
        
    D = dist2(P,Q);

    while ~isempty(notassigned)
        
        [min_v,min_i] = min(D(:));
        
        [cid, pid] = ind2sub(size(D),min_i);
        
        label(notassigned(pid)) = label(assigned(cid));
        
        P = [P;Q(pid,:)];
        Q(pid,:) = [];
        D(:,pid) = [];
       
        d = dist2(P(end,:),Q);
        D = [D;d];
        
        assigned = [assigned; notassigned(pid)];
        notassigned(pid) = [];
        
       
    end
end
