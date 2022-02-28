function out = get_subwindow(im, pos, sz)
%GET_SUBWINDOW Obtain sub-window from image, with replication-padding.
%   Returns sub-window of image IM centered at POS ([y, x] coordinates),
%   with size SZ ([height, width]). If any pixels are outside of the image,
%   they will replicate the values at the borders.
%
%   Joao F. Henriques, 2014
%   http://www.isr.uc.pt/~henriques/

if isscalar(sz)  %square sub-window方形子窗口
    sz = [sz, sz];
end
% pos 中心点坐标（y,x) (行，列）
ys = floor(pos(1)) + (1:sz(1)) - floor(sz(1)/2);%从框最左边的点，1间隔，取到最右边的点
xs = floor(pos(2)) + (1:sz(2)) - floor(sz(2)/2);%从框最上边的点，1间隔，取到最下边的点

% Check for out-of-bounds coordinates, and set them to the values at the borders
xs = clamp(xs, 1, size(im,2));
ys = clamp(ys, 1, size(im,1));

%extract image
out = im(ys, xs, :);%取了输入图像的目标框内的矩阵

end

function y = clamp(x, lb, ub)
% Clamp the value using lowerBound and upperBound

y = max(x, lb);
y = min(y, ub);

end