function window_sz = get_search_window_MCCF( target_sz, im_sz, padding)
% GET_SEARCH_WINDOW
%padding 
% if(target_sz(1)/target_sz(2) > 2)
if ((target_sz(1)/target_sz(2) > 2) && (im_sz(1)/target_sz(2) > 11))
    % For objects with large height, we restrict the search window with padding.height
    window_sz = floor(target_sz.*[1+padding.height, 1+padding.generic]);
                                    %0.4                 1.8
else
    if(prod(target_sz)/prod(im_sz(1:2)) > 0.04) %0.0662
        % For objects with large height and width and accounting for at least 10 percent of the whole image,
        % we only search 2x height and width
        window_sz=floor(target_sz*(1+padding.large));
        % 1
    else
        %otherwise, we use the padding configuration
        window_sz = floor(target_sz * (1 + padding.generic));
    end
    
end


end