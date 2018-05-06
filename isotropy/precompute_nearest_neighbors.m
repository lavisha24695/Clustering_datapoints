% Function to pre-compute K nearest neighbors

function data_details = precompute_nearest_neighbors(data, K)

%addpath vlfeat-0.9.9/toolbox/
%vl_setup();

run('C:/Users/lavis/Documents/vlfeat/vlfeat-0.9.9/toolbox/vl_setup')
vl_version verbose

fprintf('Build kd-tree\n');
data = data';
%vl_kdtreebuild;
kdtree_data = vl_kdtreebuild(data, 'verbose', 'thresholdmethod', ...
        'median');
numPts = size(data,2);
    
    for i=1:numPts

        if mod(i,100) == 0
            fprintf('.');
        end
        if mod(i,1001) == 0
            fprintf('\n');
        end

        [nbr_index, distance] = vl_kdtreequery (kdtree_data, data, ...
            data(:,i), 'numneighbors', K);
        
        data_details(i).center = data(:,i)';
        data_details(i).order = nbr_index;
        data_details(i).radius = sqrt(distance);
        
    end

    clear kdtree_data;
    fprintf('Done\n');
    
    
end