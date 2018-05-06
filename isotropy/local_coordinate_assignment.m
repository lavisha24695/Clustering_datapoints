%Function to implement local coordinate assignment in a n-dimensional
%cluster.
function [data_coordinates,data_dimensions] = local_coordinate_assignment(data,data_dimensions,data_neighbors,K,method)
%First get estimates of local coordinates using (either PCA or Roweis'
%method)
fprintf('\nLocal Coordinate Assignment\n');
numPts = size(data,1);
data_coordinates= [ ];
    for i=1:numPts
        data_coordinates(i).axes = zeros(size(data,2),data_dimensions(i));
    end
    switch method   

        case 'localPCA'
            for i=1:numPts

                if data_dimensions(i) == 0
                    continue;
                end

                if mod(i,floor(numPts/10)) == 0

                    fprintf('*');
                end

                tempD = data(data_neighbors(i).order...
                    (1:min(K*data_dimensions(i),...
                    length(data_neighbors(i).order))),:);

                cov_sample = cov(tempD);

                [Q,P]=eig(cov_sample);
                P1=P/max(diag(P));

                [vals,ord] = sort(diag(P1),'descend');


                dim = min(data_dimensions(i),size(data,2));

                for j=1:dim
                    data_coordinates(i).axes(:,j) = Q(:,ord(j)); 
                end

                Q=data_coordinates(i).axes(:,1:data_dimensions(i));

                proj_matrix = Q*Q';

                %Ortho Components is a Dxn matrix
                ortho_components = tempD*proj_matrix;
                proj_dist = tempD - ortho_components;
                proj_dist = sqrt(sum(proj_dist.*proj_dist,2));

                data_coordinates(i).residual = ...
                    mean(proj_dist)+3*std(proj_dist);

            end      
        case 'SnapShot' %SnapShot method for calculating local subspace.
            
           for i=1:numPts
                if data_dimensions(i) == 0
                    continue;
                end
                if mod(i,floor(numPts/10)) == 0
                    fprintf('*');
                end
                X=data(data_neighbors(i).order...
                    (1:min(K*data_dimensions(i),...
                    length(data_neighbors(i).order))),:);
                X=X-repmat(mean(X),size(X,1),1);

                C=X*X';
                [q,v] = eig(C);
                v1=v/max(diag(v));
                [vals,ord] = sort(diag(v1),'descend');
                dim = min(data_dimensions(i),size(data,2));
                tempQ = q(:,ord(1:dim));
                Q = X'*tempQ;
                for j=1:dim
                    data_coordinates(i).axes(:,j) = Q(:,j)/norm(Q(:,j)); 
                end
                Q=data_coordinates(i).axes(:,1:data_dimensions(i));
                proj_matrix = Q*Q';
                %Ortho Components is a Dxn matrix
                ortho_components = X*proj_matrix;
                proj_dist =X - ortho_components;
                proj_dist = sqrt(sum(proj_dist.*proj_dist,2));
                %Don't tolerate anything one standard deviation
                %away from the mean.
                %May 2009: Try different standard deviations
                data_coordinates(i).residual = ...
                    mean(proj_dist)+3*std(proj_dist);
           end
        otherwise
            fprintf(' Unknown method for subspace estimation\n\n');
            
            
    end

end
