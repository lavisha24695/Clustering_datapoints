function  [clusters, votes] = ...
    force_clustering_sum(data, data_dimensions, data_neighbors, params)

    dim =size(data,2);
    % First two columns are scalar votes (although we don't use the second
    % column any more).
    % 3:dim+2 Vector Votes (weighted)
    % dim+3:end Vector votes (uniform) (not used here)
    votes = zeros(size(data,1),2*dim+2);

    fprintf('\nComputing Shift Vectors\n');
    win_jmp  = params.Kjump;
    alpha_min = min(params.alpha_list);
    
    for i=1:size(data,1)

        temp_votes = zeros(size(data,1),1);
        diagnostics(i);

        if data_dimensions(i) == 0
            clusters.isborder_data(i) = 1;
            clusters.radius_of_uniformity(i,2)=5;
            clusters.data_neighbors(i).rnn = ...
                data_neighbors(i).order(1:clusters.radius_of_uniformity(i,2));
            clusters.data_neighbors(i).rnn_index = ...
                1:clusters.radius_of_uniformity(i,2);
            continue;
				
        end

        % Center the data
        neighbors = data(data_neighbors(i).order(1:params.Kmax),:) ...
            -repmat(data_neighbors(i).center,params.Kmax,1);
        
        win_size = params.startK;
        Kmax = params.Kmax;
        
        proj_coord = neighbors;
        clusters.radius_of_uniformity(i,2) = win_size;
        
        % Comptute force
        force = ...
            sum(proj_coord(1:win_size,:)...
            /data_neighbors(i).radius(win_size),1);
        local_dim = data_dimensions(i);
        
        % Normalized for CLT based isotropy computation
        force = force/(sqrt(win_size)*(1/sqrt(local_dim + 2)));
        
        % Compute degree of isotropy
        force_val = norm(force);
        
        % Square of the length of the force has a Chi-Squared distribution
        % with local_dim degrees of freedom.
        alpha_local = 1-chi2cdf(force_val*force_val,local_dim);
        
        temp_votes(data_neighbors(i).order(1:win_size)) = alpha_local;
            
        current_alpha_min = alpha_local;
                   
        % Repeat as long as criterion is satisfied.
        while alpha_local > alpha_min && win_size < Kmax
           
            win_size = win_size + win_jmp;
            force = ...
                sum(proj_coord(1:win_size,:)...
                /data_neighbors(i).radius(win_size),1);
            force = force/(sqrt(win_size)*(1/sqrt(local_dim + 2)));
            force_val = norm(force);
            
            alpha_local = 1-chi2cdf(force_val*force_val,local_dim);
            
            alpha_local = min(alpha_local, current_alpha_min);
            current_alpha_min = alpha_local;

            temp_votes...
                (data_neighbors(i).order(win_size-win_jmp+1:win_size)) ...
                = alpha_local;
 
        end
        
        % Compute votes
        votes_local = temp_votes(data_neighbors(i).order(1:win_size));
        
        votes(data_neighbors(i).order(1:win_size),1) = ...
            votes(data_neighbors(i).order(1:win_size),1) + votes_local;
        
        votes(data_neighbors(i).order(1:win_size),3:dim+2) = ...
            votes(data_neighbors(i).order(1:win_size),3:dim+2) ...
                        + ...
                        (repmat(votes_local,1,dim).*...
                        (repmat(data(i,:),win_size,1) ...
                        - ...
                        data(data_neighbors(i).order(1:win_size),:)))./...
                        repmat(max(eps,...
                        data_neighbors(i).radius(1:win_size)),1,dim);
                                
        clusters.radius_of_uniformity(i,2)=win_size - win_jmp;
        % rnn : neighbors over which random walks can be performed. For now
        % we use the neighborhood of isotropy. However, you can use more
        % interesting structures, e.g. relative isotropic neighbors, that
        % is neighbors that lie in each others isotropic neighborhoods.
        % Alternatively, the neighbors on which to perform the random walks
        % can be determined by the neighbors that find you in their
        % isotropic neighborhoods.
        clusters.data_neighbors(i).rnn = ...
            data_neighbors(i).order(1:clusters.radius_of_uniformity(i,2));
        clusters.data_neighbors(i).rnn_index = ...
            1:clusters.radius_of_uniformity(i,2);
        % Not used now.
        clusters.isoutlier(i) = 0;
    end
 
    fprintf('\nDone computing shift vectors\n');
        
        
end


              
                    

                        
                        
    
                        
                               
                
                        
