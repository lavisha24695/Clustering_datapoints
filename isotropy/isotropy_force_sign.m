function  pvals = ...
    isotropy_force_sign(data,...
    data_dimensions, data_coordinates, data_neighbors, params)

    pvals = zeros(size(data,1),1);
    fprintf('\nComputing PVals\n');
    
    for i=1:size(data,1)

        diagnostics(i);
        if data_dimensions(i) == 0
            continue;
        end

        neighbors = data(data_neighbors(i).order(1:params.startK),:) -...
            repmat(data_neighbors(i).center,params.startK,1);
        
        win_size = params.startK;
        win_size = min(win_size,length(data_neighbors(i).order));

        proj_coord = neighbors;
        nbd = win_size;
        Kmax = win_size;      
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
            repmat(data_neighbors(i).radius(start_idx:params.startK),...
            1,size(proj_coord,2));
        
				proj_coord = proj_coord *...
				            data_coordinates(i).axes(:,1:data_dimensions(i));

        force_val = [];

        for temp_iter = start_idx:nbd
            force_val = [force_val; ...
                norm(sum(proj_coord(start_idx:temp_iter,:),1))];
        end
        
       
        pvals(i) = signtest(...
            sign(force_val(1:end-1) - ...
            force_val(2:end)));
                    
    end
    
    fprintf('\nDone computing pvals\n');
        
        
end


              
                    

                        
                        
    
                        
                               
                
                        
