%Code to generate look-up table for the dimensionality test dims: 2-N

%function generate_cdf_look_up(total_num_dimensions,granularity,fname)

granularity = 30;
total_num_dimensions = 20;
    
    dv = 1:-2/granularity:-1;
    angles = acos(dv);
    
    syms g x;
    
    cdf=zeros(total_num_dimensions,granularity+1);
    
    for j=1:length(angles)
    
        for i=1:total_num_dimensions
        
            g=int((sin(x))^i,0,angles(j));
            
            cdf(i,j)=double(g)*gamma((i/2)+1)/(gamma((i+1)/2)*sqrt(pi));
    
        end
    
    end
 
save('look_up_cdf_31_bins.mat','cdf');    
clear g x;