function  pvals = ...
    isotropy_force_sum(data, data_dimensions, data_coordinates, data_neighbors, params)

    fprintf('\nComputing PVals\n');
    pvals = zeros(size(data,1),1);
    
    for i=1:size(data,1)

        diagnostics(i);

        if data_dimensions(i) == 0
            continue;		
        end

        % Center the data
        neighbors = data(data_neighbors(i).order(1:params.startK),:) ...
            -repmat(data_neighbors(i).center,params.startK,1);
        
        win_size = min(params.startK, size(neighbors,1));       
				proj_coord = neighbors *...
				             data_coordinates(i).axes(:,1:data_dimensions(i));

        %proj_coord = neighbors;
        
				dist = sqrt(sum(proj_coord.*proj_coord,2));

        % Comptute force
        force = ...
            sum(proj_coord(1:win_size,:)...
            /max(max(dist),eps),1);
        %force = ...
        %    sum(proj_coord(1:win_size,:)...
        %    /data_neighbors(i).radius(win_size),1);
        local_dim = data_dimensions(i);
        
        % Normalized for CLT based isotropy computation
        force = force/(sqrt(win_size)*(1/sqrt(local_dim + 2)));
        
        % Compute degree of isotropy
        force_val = norm(force);
        
        % Square of the length of the force has a Chi-Squared distribution
        % with local_dim degrees of freedom.
        pvals(i) = 1-chi2cdf(force_val*force_val,local_dim);
        
     end
 
    fprintf('\nDone pvals\n');
        
        
end


              
                    

                        
                        
    
                        
                               
                
                        
