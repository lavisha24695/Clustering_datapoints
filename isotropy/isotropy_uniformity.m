function pvals = ...
    isotropy_uniformity(data,...
    data_dimensions, ...
    data_neighbors, ...
    data_coordinates, ...
    params)
    cdf = load(params.cdf_file);
    cdf = cdf.cdf; 
    pvals = zeros(size(data,1),1);
    fprintf('\nComputing PVals\n');
    
    for i=1:size(data,1)
        diagnostics(i);
        if data_dimensions(i) == 0
            continue;
        end
        start_idx = 2;
        Kmax = params.startK;
        % Skip overlapping points
        for temp_iter = 1:Kmax   
            if data_neighbors(i).radius(temp_iter)>0
                start_idx = temp_iter;
                break;
            end 
        end   
        neighbors = data(data_neighbors(i).order(start_idx:params.startK),:) -...
            repmat(data_neighbors(i).center,params.startK-start_idx+1,1);
        
        win_size = params.startK;
        win_size = min(win_size,size(neighbors,1));

        proj_coord = neighbors *...
            data_coordinates(i).axes(:,1:data_dimensions(i));
        nbd = win_size;
        
        proj_nbd = proj_coord(1:nbd,:);
        proj_nbd = proj_nbd./...
            repmat(2*(max(data_neighbors(i).radius(nbd+start_idx-1),eps)), ...
            size(proj_nbd,1),size(proj_nbd,2)) ...
            +0.5*ones(size(proj_nbd)); 
        pval_min = 1;
        for d=1:data_dimensions(i)
           [h,pval]=kstest_local(proj_nbd(:,d),...
               [(0:1/100:1);cdf(data_dimensions(i),:)]',...
               0.01);
           if pval < pval_min
               pval_min = pval;
           end
        end
        pvals(i) = pval_min;
    end
    fprintf('\nDone computing shift vectors\n');
end
    
