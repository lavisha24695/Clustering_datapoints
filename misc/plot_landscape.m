%function to plot a landscape.

function plot_landscape(data,granularity)

%data is assumed to have (x,y,value)

X=minmax(data(:,1:2)');

x=X(1,1):(X(1,2)-X(1,1))/granularity:X(1,2);
y=X(2,1):(X(2,2)-X(2,1))/granularity:X(2,2);

Z=zeros(length(y),length(x));

for i=1:length(x)
    for j=1:length(y)
   
         eDMatrix=euclidean_distance([x(i);y(j)],data(:,1:2)');
         [eDV,I]=sort(eDMatrix,'ascend');
         Z(j,i)=mean(data(I(1:10),3));
        
    end
end

figure
surf(x,y,Z);
