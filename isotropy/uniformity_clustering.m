function [clusters, votes] = ...
    uniformity_clustering(data,...
    data_dimensions, ...
    data_neighbors, ...
    data_coordinates, ...
    params)
    
    cdf = load(params.cdf_file);
    cdf = cdf.cdf;
    
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

        start_idx = 2;
        
        Kmax = params.Kmax;
        % Skip overlapping points
        for temp_iter = 1:Kmax   
            if data_neighbors(i).radius(temp_iter)>0
                start_idx = temp_iter;
                break;
            end
            
        end
        
				neighbors = data(data_neighbors(i).order(start_idx:params.Kmax),:) -...
            repmat(data_neighbors(i).center,params.Kmax-start_idx+1,1);
        
        win_size = params.startK;
        win_size = min(win_size,length(data_neighbors(i).order));

        proj_coord = neighbors *...
            data_coordinates(i).axes(:,1:data_dimensions(i));
        nbd = win_size;
        clusters.radius_of_uniformity(i,2)= win_size;
        
        %Starting with the most stringent constraints on alpha, start
        %determining the radius of uniformity at that alpha. Once the
        %anisotropy is detected, relax alpha and find the next
        %break point.
        a=1;
        
        while nbd <= Kmax && a<=length(params.alpha_list)
            proj_nbd = proj_coord(1:nbd,:);
            proj_nbd = proj_nbd./...
                repmat(2*(max(data_neighbors(i).radius(nbd+start_idx-1),eps)), ...
                size(proj_nbd,1),size(proj_nbd,2)) ...
                +0.5*ones(size(proj_nbd)); 

            %Compare projections along each dimension with expected
            %distribution.
            isuniform=1;
                
            for d=1:data_dimensions(i)
               h=kstest_local(proj_nbd(:,d),...
                   [(0:1/100:1);cdf(data_dimensions(i),:)]',...
                   params.alpha_list(a));

                if h==1
                    isuniform=0;
                    break;
                end

            end
            
            if isuniform == 0

                vote_win = nbd + start_idx;
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
            end   
        end

        %If all ALPHA's have not been evaluated dump here.
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
    
