function [elongation, area2, perimeter2] = compute_voronoiprops(data)
    
    n = size(data,1);
    d = size(data,2);
    elongation = zeros(n,1);
    area = zeros(n,1);
    perimeter = zeros(n,1);
    [V,C] = voronoin(data);
    for i = 1:1:n
        polygonpts = V(C{i},:);
        polygonpts = polygonpts(sum(polygonpts,d)~=Inf, :);
        %%Computation of elongation
        [K, vol] = convhull(polygonpts);
        n_pts = size(polygonpts, 1);
        perim = 0;
        for j = 1:1:n_pts-1
           len = polygonpts(j, :) - polygonpts(j+1, :); 
           perim = perim + sqrt(sumsqr(len));
        end
        %length of the joining edge from the last point to the first point
        len = polygonpts(n_pts, :) - polygonpts(1, :); 
        perim = perim + sqrt(sumsqr(len));
        elongation(i) = vol/(perim*perim);
        area(i) = vol;
        area2 = area;
        mean1 = mean(area);
        std1 = std(area);
        area2(find(area>mean1 + 2*std1)) = mean1 + 2*std1;
        area2 = log(area2);
        
        perimeter(i) = perim;
        perimeter2 = perimeter;
        mean1 = mean(perimeter);
        std1 = std(perimeter);
        perimeter2(find(perimeter>mean1 + 2*std1)) = mean1 + 2*std1;
        %if(perimeter ==0)
        perimeter2 = log(perimeter2);
     end
end