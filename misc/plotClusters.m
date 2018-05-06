%Code to plot all the 2D and 3D clusters as determined by the clustering
%program. 

%Sanketh Shetty January 2007.

%data = clustered data.
%dim = dimensionality of the data.

function plotClusters(data,dim,flooze,ax2)

%sort clustered data by clusterID.
%figure;
%data;

	d=data(find(data(:,dim+1)>-1),:);
	e=data(find(data(:,dim+1)==-1),:);

	d=sortrows(d,dim+1);

	cluster_labels=unique(d(:,dim+1));

	colors='kcyrbgm';
	marks='.*x+o><*sdv^ph';

	k=1;

	G=d(:,1:dim);

	M=minmax(G');

	ax2 = [M(1,:) M(2,:)];

	ax3 = [ax2 M(dim,:)];

	mod_c=[];

	for k=1:length(cluster_labels)
	    
	    L=find(d(:,dim+1)==cluster_labels(k));
	    if length(L) > 0
	    mod_c =[mod_c;cluster_labels(k) length(L)];

	    end
	end


	if size(mod_c,1)<=100
	    LL=size(mod_c,1);
	else
	    LL=100;
	end

if length(mod_c)~=0
    
cluster_labels = sortrows(mod_c,-2);
cluster_labels(1:LL,:)

if ~flooze
    return;
end

k=1;

figure;

hold on;
if dim == 2
    plot(e(:,1),e(:,2),'k.')
end
if dim == 3
    plot3(e(:,1),e(:,2),e(:,3),'k.');
end
            
xlabel('feature 1');
ylabel('feature 2');
grid on;
%axis([-15 15 -15 15]);
while k <= length(cluster_labels)

    
    %axis([-.6 .6 -.6 .6 -.6 .6]);
    %axis(ax3);
    hold on;
    axis(ax2);
    
    for i=1:length(marks)
        for j=1:length(colors)
        
%             if colors(j) == 'k' && marks(i) == '.'
%                 continue;
%             end
            
            legend=[colors(j) marks(i)];
        
            list=find(d(:,dim+1)==cluster_labels(k));

            if dim == 2 && length(list) > 0
                
            plot(d(list,1),d(list,2),legend,'MarkerSize',20);
            k=k+1;
            end
            
            if dim == 3
            plot3(d(list,1),d(list,2),d(list,3),legend);
            k=k+1;
            end
            
           dlist=d(list,:);
             
%            dlist=sortrows(dlist,-4);
%           
            if flooze==0
                for b=1:size(dlist,1)

%                 if dlist(b,5) == 1
%                 clr = 'b';
%                 else
%                 clr = 'r';
%                 end
                
                    drawcircle(dlist(b,2),dlist(b,1),2*dlist(b,5),'r');
                %plot(dlist(b,1),dlist(b,2),'r+');
                %fprintf('%f %f %f\n',dlist(b,1),dlist(b,2),dlist(b,5));
    
                end
            end
            if k > length(cluster_labels)
                break;
            end
        end

        
            if k > length(cluster_labels)
                break;
            end

    end
    
    break;
    
end
end
%figure;
%plot(e(:,1),e(:,2),'bh');
