%% Function to estimate dimension
function data_dimension = ...
    estimate_dimensionality(data,data_neighbors,params)

method_dim = params.method_dim;
K = 30;
fprintf('Start Dimension Estimation by %s\n',method_dim);

data_dimension = zeros(size(data,1),1);

    switch method_dim

        case 'proj_L2'
            %Compare the projection along the first principal direction
            %with the expected CDF, in the L2 sense
            cdf_file = load('look_up_cdf_101_bins_50dim.mat');
            cdf = cdf_file.cdf;
            
            for i=1:size(data,1)
            
                if mod(i,100) == 0
                    fprintf('.');
                end
                if mod(i,1001) == 0
                    fprintf('\n');
                end
                neighbors = data(data_neighbors(i).order(1:K),:)-...
                    repmat(data_neighbors(i).center,K,1);
                %Obtain 1-dimensional projections
                proj = neighbors*data_neighbors(i).principal';
            
                %Obtain CDF of these projections
                proj = proj/(2*data_neighbors(i).radius(K))+0.5;
                proj = sort(proj);
                
                if proj(1) >= eps
                    
                    proj = [0;proj];
                
                end
                
                if proj(end) <= 1-eps
                
                    proj = [proj;1];
                
                end
                
                proj = unique(proj);
                
                pts_eval = [0:1/100:1];
  
                if params.smooth == 1
                    cdf_smooth = ...
                        ksdensity(proj,pts_eval,'width',1,'function','cdf');
                    cdf_dist = ...
                        sum((cdf - repmat(cdf_smooth,size(cdf,1),1)).^2,2);
                    [dummy,data_dimension(i)] = min(cdf_dist);
                    
                else
                     cdf_regular = 0:1/(size(proj,1)-1):1;                         
                     cdf_regular = interp1(proj,...
                         cdf_regular,pts_eval,'linear');
                     cdf_dist = ...
                         sum((cdf - ...
                         repmat(cdf_regular,size(cdf,1),1)).^2,2);
                     [dummy,data_dimension(i)] = min(cdf_dist);
                    
                end
            
            end
            
            
        case 'proj_L1'
            
            %Compare the projection along the first principal direction
            %with the expected CDF, in the L1 sense 
            cdf_file = load('look_up_cdf_101_bins_50dim.mat');
            cdf = cdf_file.cdf;
            
            for i=1:size(data,1)
            
                
                if mod(i,100) == 0
                    fprintf('.');
                end
                if mod(i,1001) == 0
                    fprintf('\n');
                end
                neighbors = data(data_neighbors(i).order,:)...
                    -repmat(data_neighbors(i).center,K,1);
                %Obtain 1-dimensional projections
                proj = neighbors*data_neighbors(i).principal';
            
                %Obtain CDF of these projections
                proj = proj/(2*data_neighbors(i).radius(K))+0.5;
                proj = sort(proj);
                
                if proj(1) >= eps
                    
                    proj = [0;proj];
                
                end
                
                if proj(end) <= 1-eps
                
                    proj = [proj;1];
                
                end
                
                proj = unique(proj);
                
                pts_eval = [0:1/100:1];
                   
                if params.smooth == 1
                    cdf_smooth = ...
                        ksdensity...
                        (proj,pts_eval,'width',1,'function','cdf');
                    cdf_dist = ...
                        sum(abs(cdf - repmat(cdf_smooth,size(cdf,1),1)),2);
                    [dummy,data_dimension(i)] = min(cdf_dist);
                    
                else

                    cdf_regular = 0:1/(size(proj,1)-1):1;
                    cdf_regular = ...
                        interp1(proj,cdf_regular,pts_eval,'linear');
                    cdf_dist = ...
                        sum(abs(cdf - repmat(cdf_regular,size(cdf,1),1)),2);
                    [dummy,data_dimension(i)] = min(cdf_dist);
                    
                end
            end
          
       
        case 'MLE'
       
            % Maximum Likelihood Estimate of local dimensions based on 
            % Poisson Data distribution assumption 
            k1 = min(10,size(data,1)-10);
            k2 = min(20,size(data,1)-2);
            
            for i=1:size(data,1)
            
                
                if mod(i,100) == 0
                    fprintf('.');
                end
                if mod(i,1001) == 0
                    fprintf('\n');
                end
                
                neighbors = data(data_neighbors(i).order(2:k2+1),:);
                                
                dist = euclidean_distance(data(i,:)',neighbors');     
                
                temp = log(dist);
                temp2 = cumsum(temp);
                
                %October 19 2008: Fixed error in MLE dimension computation
                num1 = temp(k1:k2);
                num2 = temp2(k1-1:k2-1)./(k1-1:k2-1);
                
                tot_num = 1./(num1-num2);
                data_dimension(i) = sum(tot_num)/(k2-k1+1);
            end
            
            
        case 'NearestNbr_dim'
        
            % Set neighborhood range to search in
            k1 = 6;
            k2 = 12;
            
            % Compute nearest neighbor based dimension 
            
            for i=1:size(data,1)
                
                
                if mod(i,100) == 0
                    fprintf('.');
                end
                if mod(i,1001) == 0
                    fprintf('\n');
                end
            
                neighbors = data(data_neighbors(i).order(1:k2+1),:);
                dist = euclidean_distance(data(i,:)',neighbors');     
                data_dimension(i) = ...
                    (log(dist(k2+1))-log(dist(k1+1)))/(log(k2)-log(k1));
           
            end
            
        case 'Eigen'
        
            %Considering this is a local eigen value problem we need to the
            % quality of the local PCA goes down with decreasing number of
            % points.
            for i=1:size(data,1)
                
                
                if mod(i,100) == 0
                    fprintf('.');
                end
                if mod(i,1001) == 0
                    fprintf('\n');
                end
                if size(data,2) < 30
                    neighbors = data(data_neighbors(i).order,:);   
                    cov_nbr = cov(neighbors);
                else
                    eDMatrix=euclidean_distance(data(i,:)',data'); 
                    [eDV,I]=sort(eDMatrix,'ascend');
                    cov_nbr = cov(data(I(1:size(data,2)*10)));
                end

                [Q,P] = eig(cov_nbr);
                P = P/max(diag(P));
                
                data_dimension(i) = sum(diag(P)>0.3);
                
            end
            
            
        otherwise

            fprintf('Unknown dimension estimator!!!\n.Quitting.\n');
    end

fprintf('\n');


end
