function  [clusters, votes] = ...
    force_clustering_sign(data,...
    data_dimensions, data_neighbors, params)

    dim = size(data,2);
    votes = zeros(size(data,1),2*dim+2);

    fprintf('\nComputing Shift Vectors\n');
    
    for i=1:size(data,1)

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

        neighbors = data(data_neighbors(i).order(1:params.Kmax),:) -...
            repmat(data_neighbors(i).center,params.Kmax,1);
        
        win_size = params.startK;
        win_size = min(win_size,length(data_neighbors(i).order));

        Kmax = params.Kmax;
        
        proj_coord = neighbors;
        nbd = win_size;
        clusters.radius_of_uniformity(i,2)= win_size;
        
        %Starting with the most stringent constraints on alpha, start
        %determining the radius of uniformity at that alpha. Once the
        %anisotropy is detected, relax alpha and find the next
        %break point.
        a=1;
        
        start_idx = 2;
        
        % Skip overlapping points
        for temp_iter = 1:Kmax   
            if data_neighbors(i).radius(temp_iter)>0
                start_idx = temp_iter;
                break;
            end
            
        end
        
        proj_coord(start_idx:end,:) = ...
            proj_coord(start_idx:end,:)./...
            repmat(data_neighbors(i).radius(start_idx:end),...
            1,size(proj_coord,2));

        force_val = [];

        for temp_iter = start_idx:nbd+start_idx-1
            force_val = [force_val; ...
                norm(sum(proj_coord(start_idx:temp_iter,:),1))];
        end
        
        force = sum(proj_coord(start_idx:nbd+start_idx-1,:),1);
        
        while nbd+start_idx <= Kmax && a<=length(params.alpha_list)
                
            % Check isotropy
            if length(force_val)<=win_size
                isisotropic = 1;
            else

                pval = signtest(...
                    sign(force_val(end-win_size:end-1) - ...
                    force_val(end-win_size+1:end)));
                isisotropic = pval>params.alpha_list(a);
                    
            end
            
            % Vote if isotropy criterion is violated
            if isisotropic == 0

                vote_win = nbd;
                if a==1
                    clusters.radius_of_uniformity(i,2)=vote_win;
                end
                % Do not vote till at least one neighborhood has been found 
                % uniform at some alpha level. 
                % If none of the alpha's pan out vote once at the 
                % neighborhood size win_size.	
                if nbd == win_size && a<length(params.alpha_list)
                    a=a+1;
                    continue;	
                end

                clusters.radius_of_uniformity(i,1)...
                    =data_neighbors(i).radius(vote_win);
                
                %Scalar votes
                wts=(1-(data_neighbors(i).radius(1:vote_win)...
                    /data_neighbors(i).radius(vote_win)));
                votes(data_neighbors(i).order(1:vote_win),2) = ...
                    votes(data_neighbors(i).order(1:vote_win),2)+1;
                votes(data_neighbors(i).order(1:vote_win),1) = ...
                    votes(data_neighbors(i).order(1:vote_win),1)+wts;

                %Vector votes
                % Note: Bug fixed. The votes were not normalized earlier.
                votes(data_neighbors(i).order(1:vote_win),3:dim+2) = ...
                    votes(data_neighbors(i).order(1:vote_win),3:dim+2) ...
                    +(repmat(wts,1,dim).*(repmat(data(i,:),vote_win,1)...
                    - data(data_neighbors(i).order(1:vote_win),:)))./...
                    repmat(...
                    max(eps,data_neighbors(i).radius(1:vote_win)),...
                    1,dim);

                votes(data_neighbors(i).order(1:vote_win),...
                    dim+3:2*dim+2) = ...
                    votes(data_neighbors(i).order(1:vote_win), ...
                    dim+3:2*dim+2) ...
                    +(repmat(data(i,:),vote_win,1) ...
                    -data(data_neighbors(i).order(1:vote_win),:))./...
                    repmat(max(eps,...
                    data_neighbors(i).radius(1:vote_win)),1,dim);

                %Do not increment ring value. Instead increment the index
                %of the significance value.
                a=a+1;
            else
                nbd=nbd+params.Kjump;
                % Compute force over new neighborhood
                force = force +...
                    sum(proj_coord...
                    (nbd-params.Kjump+start_idx:nbd+start_idx,:),1);
                force_val = [force_val; norm(force)];

            end   
        end

        %If all alpha's have not been evaluated dump here.
        if a < length(params.alpha_list)

            factor = length(params.alpha_list)-a;

            if nbd~=win_size
                vote_win = min(nbd,length(data_neighbors(i).order));
            else
                vote_win = min(nbd,length(data_neighbors(i).order));
            end

            clusters.radius_of_uniformity(i,1) =...
                data_neighbors(i).radius(vote_win);
            if a==1
                clusters.radius_of_uniformity(i,2)=vote_win;
            end

            %Scalar votes
            wts=factor*(1-(data_neighbors(i).radius(1:vote_win)...
                /data_neighbors(i).radius(vote_win)));
            votes(data_neighbors(i).order(1:vote_win),2) = ...
                votes(data_neighbors(i).order(1:vote_win),2)+1;
            votes(data_neighbors(i).order(1:vote_win),1) = ...
                votes(data_neighbors(i).order(1:vote_win),1)+wts;

            
            %Vector votes
            votes(data_neighbors(i).order(1:vote_win),2) = ...
                    votes(data_neighbors(i).order(1:vote_win),2)+1;
                votes(data_neighbors(i).order(1:vote_win),1) = ...
                    votes(data_neighbors(i).order(1:vote_win),1)+wts;

            %Vector votes
            % Note: Bug fixed. The votes were not normalized earlier.
            votes(data_neighbors(i).order(1:vote_win),3:dim+2) = ...
                votes(data_neighbors(i).order(1:vote_win),3:dim+2) ...
                +(repmat(wts,1,dim).*(repmat(data(i,:),vote_win,1)...
                - data(data_neighbors(i).order(1:vote_win),:)))./...
                repmat(...
                max(eps,data_neighbors(i).radius(1:vote_win)),...
                1,dim);

            votes(data_neighbors(i).order(1:vote_win),...
                dim+3:2*dim+2) = ...
                votes(data_neighbors(i).order(1:vote_win), ...
                dim+3:2*dim+2) ...
                +(repmat(data(i,:),vote_win,1) ...
                -data(data_neighbors(i).order(1:vote_win),:))./...
                repmat(max(eps,...
                data_neighbors(i).radius(1:vote_win)),1,dim);
        end
        
        clusters.radius_of_uniformity(i,2)= nbd + start_idx - params.Kjump;
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


              
                    

                        
                        
    
                        
                               
                
                        
