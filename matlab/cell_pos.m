

% find the positions of all cells in the mask
% kieran.campbell@sjc.ox.ac.uk

function[xy] = cell_pos(mask_xell)

mask_neg = ~(mask_xell);
CC = bwconncomp(mask_neg);


N = length(CC.PixelIdxList);

xy = zeros(N,2);

for i = 1:N,
   mask_neg_1 = zeros(size(mask_neg));
   mask_neg_1(CC.PixelIdxList{i}) = 1;
   [x, y] = center_of_mass(mask_neg_1);
   xy(i, 1) = x;
   xy(i, 2) = y;
end

end