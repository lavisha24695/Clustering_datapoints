% Author : Sanketh Shetty
function [pvals, data_neighbors, data_dimensions] = ...
    isotropy_clustering(data, criterion, win_size)

    % Parameters:
    % Change these for speedups. They do not have much of an effect on
    % clustering.
    params.Kmax = min(200, size(data,1)); % maximum window size for testing
    % Subspace estimation
    % options: 'MLE','Eigen','NearestNbr_dim','proj_L1','proj_L2',
    params.method_dim = 'MLE';
    %options: 'SnapShot','localPCA';
    params.subspace = 'SnapShot';
    params.subspace_refine = 0; % 1: Refines neighborhood points to include 
    % only those in the linear subspace
    params.Kper_dim = 10; % Number of points x dimension gives number of points
    % used in estimating the local linear subspace
    params.Krefine_dim = 2;
    params.method_refine = 'mode';
    params.Kper_dimC = 10;
   	params.cdf_file = 'isotropy/look_up_cdf_101_bins_50dim.mat';
    params.startK = min(win_size, size(data,1));

    if strcmp(criterion, 'force-sum')
        
        % Precompute Neighbors
        data_neighbors = precompute_nearest_neighbors(data, params.Kmax);
        
        %Estimate intrinsic dimensions of point neighborhoods
        data_dimensions = ...
            estimate_dimensionality(data,data_neighbors,params);

        %Refine the dimension estimates over local neighborhoods
        data_dimensions = ...
            refine_dimension_estimates(data_dimensions, data_neighbors,...
            params.method_refine,params.Krefine_dim);
        data_dimensions = min(round(data_dimensions), size(data,2));
       
        %Assign Local Coordinates
        [data_coordinates,data_dimensions] = ...
            local_coordinate_assignment(data,data_dimensions,...
            data_neighbors,params.Kper_dim,params.subspace);        
        
				pvals = isotropy_force_sum(data, ...
            data_dimensions, data_coordinates, data_neighbors, params);
    end

    if strcmp(criterion, 'force-sign')
        % Precompute Neighbors
        data_neighbors = precompute_nearest_neighbors(data, params.Kmax);
        
        %Estimate intrinsic dimensions of point neighborhoods
        data_dimensions = ...
            estimate_dimensionality...
            (data,data_neighbors,params);

         %Refine the dimension estimates over local neighborhoods
         data_dimensions = ...
            refine_dimension_estimates(data_dimensions, data_neighbors,...
            params.method_refine,params.Krefine_dim);
         data_dimensions = min(round(data_dimensions), size(data,2));
        
        %Assign Local Coordinates
        [data_coordinates,data_dimensions] = ...
            local_coordinate_assignment(data,data_dimensions,...
            data_neighbors,params.Kper_dim,params.subspace);        
         
				 pvals = isotropy_force_sign(data, ...
            data_dimensions, data_coordinates, data_neighbors, params);
    end

    if strcmp(criterion, 'uniformity')
        % Precompute Neighbors
        data_neighbors = precompute_nearest_neighbors(data, params.Kmax);        
        
        %Estimate intrinsic dimensions of point neighborhoods
        data_dimensions = ...
            estimate_dimensionality...
            (data,data_neighbors,params);

        %Refine the dimension estimates over local neighborhoods
        data_dimensions = ...
            refine_dimension_estimates(data_dimensions, data_neighbors,...
            params.method_refine,params.Krefine_dim);
        
				data_dimensions = min(round(data_dimensions), 50);
				data_dimensions = min(round(data_dimensions), size(data,2));


        %Assign Local Coordinates
        [data_coordinates,data_dimensions] = ...
            local_coordinate_assignment(data,data_dimensions,...
            data_neighbors,params.Kper_dim,params.subspace);        
        
        pvals = isotropy_uniformity(data, ...
            data_dimensions, data_neighbors, data_coordinates, params);
        
    end
end
