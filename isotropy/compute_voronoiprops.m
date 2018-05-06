function [elongation, area, perimeter] = compute_voronoiprops(data)
    
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
        perimeter(i) = perim;
     end
end