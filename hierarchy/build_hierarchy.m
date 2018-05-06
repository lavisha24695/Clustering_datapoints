%This piece of code builds a hierarchy after the inital partitioning is
%done. It first identifies merge points and then finds the persistence
%level for these merge points. The merge points are sorted in increasing
%order of persistence and then a hierarchy is generated.

function [label_evolution,plevels] = build_hierarchy(initial_labels, ...
    data_neighbors, data_dimensions, neighborhood_size)

fprintf('Building Hierarchy\n');
%Estimate density at all the points.

num_pts = length(initial_labels);
density = zeros(num_pts,1);

for i = 1:num_pts
	num_nbrs = min(neighborhood_size,length(data_neighbors(i).order));
	nbrs = data_neighbors(i).order(1:num_nbrs);
	dists = data_neighbors(i).radius(1:num_nbrs);
	dists = dists(initial_labels(nbrs)==initial_labels(i));
	nbrs = nbrs(initial_labels(nbrs)==initial_labels(i));
	
	num_nbrs = length(nbrs);	
	t_k = max(dists);
	
	if t_k == 0 || num_nbrs == 0
		density(i) = 0;	
	else
		density(i) = num_nbrs*gamma(data_dimensions(i)+0.5)/...
            (pi^(data_dimensions(i)/2)*t_k^(data_dimensions(i)));	
	end
end

unique_labels = unique(initial_labels);
density_label = zeros(length(unique_labels),1);
cluster_count = zeros(length(unique_labels),1);

for n = 1:length(unique_labels)
	selector = initial_labels == unique_labels(n);
	density_label(n) = max(density(selector));
	cluster_count(n) = sum(selector);
end

% Uncomment this and comment the next two lines following this if you have
% a lot of atomic clusters
%boundary_density = sparse(length(unique_labels),length(unique_labels));
%bdr_count = sparse(length(unique_labels),length(unique_labels));
boundary_density = zeros(length(unique_labels),length(unique_labels));
bdr_count = zeros(length(unique_labels),length(unique_labels));


for i = 1:num_pts
	num_nbrs = min(neighborhood_size,length(data_neighbors(i).order));
	index = data_neighbors(i).order(1:num_nbrs);	
	labels_in_nbd = initial_labels(index);
	points = index(labels_in_nbd~=initial_labels(i));
    if numel(points) == 0
        continue;
    end
	labels_in_nbd = labels_in_nbd(labels_in_nbd~=initial_labels(i));	
	
    for j =1:length(points)
        merge_point = points(j);
        num_nbrs2 = min(neighborhood_size,...
            length(data_neighbors(merge_point).order));
        index2 = data_neighbors(merge_point).order(1:num_nbrs2);
        
        if sum(index2 == i) == 0
            continue;
        end

        index_merge_label = find(unique_labels == labels_in_nbd(j));
        index_org = find(unique_labels == initial_labels(i));

        boundary_density(index_org,index_merge_label) = ...
            boundary_density(index_org,index_merge_label) ...
            + min(density(i),density(merge_point));
        bdr_count(index_org,index_merge_label) = ...
            bdr_count(index_org,index_merge_label) + 1;

    end
end
boundary_density = boundary_density + boundary_density';
bdr_count = bdr_count + bdr_count';

label_evolution{1} = initial_labels;
level = 1;
labels = initial_labels;

hierarchy = initialize_hierarchy(boundary_density, ...
    bdr_count, density_label, unique_labels);

hierarchy = sortrows(hierarchy,7);
next_label = max(initial_labels) + 1;

plevels = [];

    while ~isempty(hierarchy)
        next_merge = hierarchy(1,:);
        hierarchy(1,:) = [];
        change_labels = find(labels==next_merge(1)); 
        labels(change_labels) = next_label;
        change_labels = find(labels==next_merge(2));
        labels(change_labels) = next_label;

        max_den = max(next_merge(3),next_merge(4));

        plevels = [plevels;next_merge(end)];

        hierarchy = ...
                update_hierarchy_map(hierarchy,next_merge(1), ...
                next_merge(2),max_den,next_label);

        next_label = next_label + 1;
        level = level + 1;
        label_evolution{level} = labels;
    end

fprintf('Done\n');

end

function hierarchy_map = initialize_hierarchy(boundary_density, ...
    bdr_count, density_label, labels)

    densities = density_label;    

    NR = length(density_label);
    hierarchy_map = [];
    for i = 1: NR-1
        row = bdr_count(i,i+1:end);
        brow = boundary_density(i,i+1:end);
        
        if isempty(row)
            continue;
        end
    
        nbrs = find(row>0);
        
        if isempty(nbrs)
            continue;
        end
        
        den_o = repmat(densities(i),length(nbrs),1);

        min_den = max(0,min(den_o,densities(i+nbrs)));
        persistence = ...
            max(0,(min_den - (full(brow(nbrs))'./full(row(nbrs))')))./...
            max(eps,min_den);

        hierarchy_map = [hierarchy_map; ...
            repmat(labels(i),length(nbrs),1) labels(i+nbrs) ...
            repmat(density_label(i),length(nbrs),1) density_label(i+nbrs) ...
            full(row(nbrs))' full(brow(nbrs))' persistence];
    end
end

function updated_hierarchy = ...
    update_hierarchy_map(hierarchy_map,label1,label2,max_den,next_label)
    
    [c1,row1] = find_candidates(hierarchy_map,label1);
    [c2,row2] = find_candidates(hierarchy_map,label2);
    
    candidates = [c1;c2];
    hierarchy_map([row1;row2],:) = [];
    candidates = sortrows(candidates,1);
    
    i = 1;
    
    while ~isempty(candidates)
        if size(candidates,1) > 1 && candidates(i,1) == candidates(i+1,1)
            min_e = min(max_den,candidates(i,2));
            valley = sum(candidates(i:i+1,4))/sum(candidates(i:i+1,3));
            
            persistence = max(0,(min_e - valley))/max(eps,min_e);            

            hierarchy_map = [hierarchy_map; ...
                next_label candidates(i,1) max_den candidates(i,2) ...
                sum(candidates(i:i+1,3)) sum(candidates(i:i+1,4)) persistence];
            candidates(i:i+1,:) = [];
        else
            
            min_e = min(max_den,candidates(i,2));
            persistence = max(0,(min_e - candidates(i,4)/candidates(i,3)))/max(eps,min_e);
            hierarchy_map = [hierarchy_map;next_label candidates(i,1) ...
                max_den candidates(i,2) candidates(i,3) candidates(i,4) ...
                persistence];
            candidates(i,:) = [];
        end
    end
    updated_hierarchy = sortrows(hierarchy_map,7);
end


function [candidates,rows] = find_candidates(hierarchy_map,label)
    
    r2 = find(hierarchy_map(:,2) == label);
    candidates = hierarchy_map(r2,[1 3 5 6]);
    r1 = find(hierarchy_map(:,1) == label);
    % Reorder columns so that the entires correspond
    candidates = [candidates; hierarchy_map(r1,[2 4 5 6])];
    rows = [r1(:);r2(:)];
    
end



