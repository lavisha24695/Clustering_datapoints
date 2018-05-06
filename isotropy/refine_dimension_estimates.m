% Function to refine dimension estimates

function data_dimensions_refined = ...
    refine_dimension_estimates...
    (data_dimensions,data_neighbors,method_refine,K)

data_dimensions_refined = zeros(size(data_dimensions));

fprintf('Refine Dimension Estimates\n');

    switch method_refine

            case 'mode'

                for i=1:length(data_dimensions)

                    if mod(i,100) == 0
                        fprintf('.');
                    end
                    if mod(i,1001) == 0
                        fprintf('\n');
                    end
                
                    data_dimensions_refined(i) = ...
                        mode(data_dimensions(data_neighbors(i).order(1:K)));
                    
                end


            case 'median'

                
                for i=1:length(data_dimensions)

                    if mod(i,100) == 0
                        fprintf('.');
                    end
                    if mod(i,1001) == 0
                        fprintf('\n');
                    end
                
                    data_dimensions_refined(i) =...
                        median(data_dimensions(data_neighbors(i).order(1:K)));
                end


            case 'mean'

                for i=1:length(data_dimensions)

                    if mod(i,100) == 0
                        fprintf('.');
                    end
                    if mod(i,1001) == 0
                        fprintf('\n');
                    end
                
                    data_dimensions_refined(i) =...
                        mean(data_dimensions(data_neighbors(i).order(1:K)));
                end


            otherwise

                fprintf('Unknown dimension estimator!!!\n.Quitting.\n');

    end



end