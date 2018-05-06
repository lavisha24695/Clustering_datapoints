function data_votes_out = refine_data_votes(data_votes,data_coordinates,dim)

%Function that projects the directional votes in the local subspace.

fprintf('Refining Directional Votes!\n');

data_votes_out = zeros(size(data_votes));

for i=1:size(data_votes,1)

    diagnostics(i);
    %Copy vote values over.
    data_votes_out(i,1:2) = data_votes(i,1:2);
    
    data_votes_out(i,3:dim+2) = ((data_coordinates(i).axes*data_coordinates(i).axes')*data_votes(i,3:dim+2)')';
    %data_votes_out(i,3:dim+2) = ((eye(dim)-data_coordinates(i).axes*data_coordinates(i).axes')*data_votes(i,3:dim+2)')';
    data_votes_out(i,dim+3:2*dim+2) = ((data_coordinates(i).axes*data_coordinates(i).axes')*data_votes(i,dim+3:2*dim+2)')';
    %data_votes_out(i,dim+3:2*dim+2) = ((eye(dim)-data_coordinates(i).axes*data_coordinates(i).axes')*data_votes(i,dim+3:2*dim+2)')';
end

fprintf('Done Refining Directional Votes\n');


end
