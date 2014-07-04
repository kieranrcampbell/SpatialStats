% find boundary weights:
% for each cell returns a vector of length (number of nearest neighbours)
% where each entry is the number of shared pixels along the boundary
% between the two cells
% kieranrcampbell@gmail.com



function[RelBoundarySize] = get_boundary_size(PixelBoundaryList, Xell_nearest)
RelBoundarySize = cell(size(PixelBoundaryList));

for i = 1:length(PixelBoundaryList)
    nearest_cells = Xell_nearest{i}(:,1);
    for j = nearest_cells'

        cell_intersect = length(intersect(PixelBoundaryList{i}, PixelBoundaryList{j}));
        
        if(cell_intersect == 0)
            error('Nearest neighbour finding gone wrong')
        end
        RelBoundarySize{i} = [RelBoundarySize{i} cell_intersect];
    end
end

end



