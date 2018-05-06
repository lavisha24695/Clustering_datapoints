% Obsolete code (as of ver 1.1)
function [data_neighbors, data_dimensions, data_coordinates] = ...
    dim_and_coord(data,params)

    % This parameter is for 'proj_L2' and 'proj_L1' and they are not yet
    % proven as good estimators of dimension.
    K=30;
    
    Krefine_dim = params.Krefine_dim;
    Kper_dim = params.Kper_dim;
    method_dim = params.method_dim;

    %Compute Neighbors
    data_neighbors = precompute_nearest_neighbors(data, params.Kmax);

    %Estimate intrinsic dimensions of point neighborhoods
    data_dimensions = ...
        estimate_dimensionality(data,data_neighbors,K,method_dim,params);

    %Refine the dimension estimates over local neighborhoods
    data_dimensions2 = ...
        refine_dimension_estimates(data_dimensions,data_neighbors,...
        params.method_refine,Krefine_dim);

    data_dimensions2 = ...
        round(min(data_dimensions2,repmat(size(data,2),size(data,1),1)));
    
    %Assign Local Coordinates
    [data_coordinates,data_dimensions] = ...
        local_coordinate_assignment(data,data_dimensions2,...
        data_neighbors,Kper_dim,params.subspace);

    
end