

function[weights] = get_weights(mask_xell, PixelIdxList, Xell_nearest)

PixelBoundaryList = find_boundary(mask_xell, PixelIdxList);
weights = get_boundary_size(PixelBoundaryList, Xell_nearest);

end