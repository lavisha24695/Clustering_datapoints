%This piece of code is for partitioning the dataset.
%Version 1: The gradient estimated from the neighborhood is used to decide
%search direction.

%Version 2: Two different cases are taken. In case the point has a radius
%of uniformity > 0, find all points in that neighborhood that also consider
%this point in their neighborhood of uniformity. Connect all these points
%densely. 
%In case the point does not have a radius of uniformity of its own, then
%find the nearest neighbor in the direction in which its membership vector
%points.

%Version 3: In this we first use the membership vector to find a coarse
%initial partitioning of the clusters. Then we pass the partitions to the
%hierarchy building algorithm.

%Version 4: For ICCV 2009. For each point we have a preselected set of
%neighbors. Among these neighbors connect the point to the nearest neighbor in the direction
%of increasing "gradient". 

function [Labels,available_ID,border] = partition_data_set_use_membership(data_with_votes,radius_of_uniformity,data_neighbors,isoutlier,fname)

available_ID=0;

data = load(fname);
numPts = size(data,1);
cd_Labels = -ones(size(data,1),1);
border = ones(numPts,1);
parent = -ones(numPts,1);
dim = size(data,2);

evaluated = zeros(size(data,1),1);


%Label points with zero weights first. They can be outliers/cluster
%centers/valley points.
k=find(sum(data_with_votes(:,3:dim+2)==0,2)==dim);
cd_Labels(k) = 1:length(k);
assignLabel = length(k)+1;
evaluated(k) = 1;


%While there remain points to be evaluated
while length(evaluated(evaluated==0)) > 0

    unvisited = find(evaluated==0);
    i=unvisited(1);

    %Book keeping
    if cd_Labels(i) == -1 
        cd_Labels(i)=assignLabel;
        assignLabel = assignLabel+1;
    end

    if isoutlier(i)
    
        change_parent = find(parent==i);
        evaluated(change_parent) = 0;
        parent(change_parent) = -1;
        parent(i)=-2;
        evaluated(i) = 1;
        continue;
    end
    
    visited_list = i;
    reached_end = 0;
    bounce = 0;
    
    %Climb the gradient till the top of "density hill" is reached.
    while ~reached_end 

            %Mark current data point as evaluated.
            evaluated(i)=1;
            unvisited = find(evaluated~=1);

            %if mod(length(unvisited),10) == 0
            %    fprintf('\nNodes left:%d',length(unvisited));
            %end

            %data_neighbors.rnn has the relevant nearest neighbors
            %Among these relevant nearest neighbors we need to find ones in the
            %direction of the increasing "force".

            nbd_incl = data_neighbors(i).rnn(2:end);
            nbd_index = data_neighbors(i).rnn_index(2:end);
            
            no_uniform = 0;

            if isempty(nbd_incl)

                %Mark as an outlier and abandon
                reached_end = 1;
                continue;
            end
         
            %Set trigger to fail.
            reached_end = 1;
            %Calculate norm of the vector weight
            mag_den = norm(data_with_votes(i,3:dim+2));
            mag_nbd = sqrt(sum(data_with_votes(nbd_incl,3:dim+2).*data_with_votes(nbd_incl,3:dim+2),2));

            if mag_den == 0
                mag_den = eps;
            end
            
            mag_nbd(mag_nbd==0)=eps;
            candidates = data(nbd_incl,:);
            
            %We want to find points in the positive half-space defined by
            %the "gradient vector" at the point. 
            angle_vec_grad = (candidates(:,:)-repmat(data(i,:),size(candidates,1),1)).*repmat(data_with_votes(i,3:dim+2),size(candidates,1),1);
            angle_vec_grad = (sum(angle_vec_grad,2)./data_neighbors(i).radius(nbd_index))/mag_den;
            
            angle_grad_grad = data_with_votes(nbd_incl,3:dim+2).*repmat(data_with_votes(i,3:dim+2),size(candidates,1),1);
            angle_grad_grad = (sum(angle_grad_grad,2)./mag_nbd)/mag_den;
         
            %If there are no points in this direction label as outlier and
            %initiate alternate paths
             if isempty(angle_vec_grad(angle_vec_grad>0))
                 change_parent = find(parent==i);
                 evaluated(change_parent) = 0;
                 parent(change_parent) = -1;
                 parent(i)=-2;
                 
                 %evaluated(i)=0;
                 break;
             end
         
            % Determine the average behavior of the neighborhood. 
            angle_avg = mean(angle_grad_grad(angle_vec_grad>0));
          
  
            %For each candidate
             for j=1:size(candidates,1)


                 %If the point lies in the positive halfspace & hasn't
                 %previously been visited & isn't an outlier
                  if angle_vec_grad(j) > 0 && isempty(visited_list(visited_list==nbd_incl(j))) && parent(nbd_incl(j)) ~=-2 %&& data_with_votes(I(j),1)>=data_with_votes(i,1)-0.05*data_with_votes(i,1)
                  
                      %Link them
                      parent(i)=nbd_incl(j);
                      
                      %Now if this point's direction vector is pointing
                      %opposite to the current point's, it means they are
                      %converging. Move to the point, bounce and stop.
                
                      if angle_grad_grad(j) <= 0 && angle_avg <=0

                        %This keeps track of how many times the climber
                        %bounces. This happens at the peaks.
                        bounce=bounce+1;
              
                          if cd_Labels(nbd_incl(j))~=-1
                            change = find(cd_Labels==cd_Labels(nbd_incl(j)));
                            cd_Labels(change) = cd_Labels(i);
                          else
                              cd_Labels(nbd_incl(j))= cd_Labels(i);
                          end
    %                       j=j+1;
    %                       end
                           if bounce <= 1
                            reached_end = 0;
                            i=nbd_incl(j);
                            visited_list = [visited_list;nbd_incl(j)];
                           else
                               reached_end=1;
                               bounce =0;
                           end

                           break;

                      end

                      if angle_grad_grad(j) > 0 %|| angle_grad_grad(j) <=0 %&& angle_avg > 0

                           border(nbd_incl(j))=0;
                            if cd_Labels(nbd_incl(j)) ~= -1
                                reached_end = 1;
                                change = find(cd_Labels==cd_Labels(nbd_incl(j)));
                                cd_Labels(change) = cd_Labels(i);
                           
                                break;
                            else
                                reached_end = 0;
                                cd_Labels(nbd_incl(j))= cd_Labels(i);
                                visited_list = [visited_list;nbd_incl(j)];
                            
                                i=nbd_incl(j);
                                evaluated(i)=1;
                                break;
                            end    
                      end
                      
                      
                  end

              end %for
        
        if j>size(candidates,1)    
           reached_end=1;
        end
        
    end%while

end%outer while





Labels = cd_Labels;

border=border==1;

