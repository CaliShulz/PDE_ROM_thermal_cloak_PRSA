function [reduced_elements,reduced_boundaries] = get_reduced_mesh(MESH,reduced_indexes)


% Get elements in the mesh corresponding to reduced indexes and remap
% elements

% remap also boundary elements


% Eliminate elements that have nodes not included in reduced_indexes
vertices_index = 1:length(MESH.vertices);
not_reduced_basis_index = ones(1,length(MESH.vertices));
not_reduced_basis_index(reduced_indexes) = 0;
not_reduced_vertices    = vertices_index(not_reduced_basis_index == 1);
reduced_elements        = MESH.elements;

[sharedvals,index_row_1] = ismember(reduced_elements(1,:),not_reduced_vertices);
index_row_1 = find(index_row_1 ~= 0);
reduced_elements(:,index_row_1) = [];

[sharedvals,index_row_2] = ismember(reduced_elements(2,:),not_reduced_vertices);
index_row_2 = find(index_row_2 ~= 0);
reduced_elements(:,index_row_2) = [];

[sharedvals,index_row_3] = ismember(reduced_elements(3,:),not_reduced_vertices);
index_row_3 = find(index_row_3 ~= 0);
reduced_elements(:,index_row_3) = [];

reduced_boundaries = MESH.boundaries;
reduced_boundaries_small = reduced_boundaries(1:2,:);

% Remap elements in reduced mesh
% remap elements
new_index = 1;
for old_index = reduced_indexes
    
    reduced_elements(find(reduced_elements==old_index)) = new_index;
    reduced_boundaries_small(find(reduced_boundaries_small==old_index)) = new_index;
    
    new_index = new_index+1;
    
end

reduced_boundaries(1:2,:) = reduced_boundaries_small;



end

