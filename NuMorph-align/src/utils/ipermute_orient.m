function img_permuted = ipermute_orient(img, orient)
    % Inverse permute and flip to desired orientation
    % Inputs img: image array
    %        orient: row vector of dimension order, negative values
    %        indicate flip of axis
    flip_orient = orient([2 1 3]);
    disp(orient);
    flip_dim = find(flip_orient<0);
    if ~isempty(flip_dim)
        for i = flip_dim
            img= flip(img, i);
        end
    end
    img_permuted = ipermute(img, abs(orient));
end