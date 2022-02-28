function labels = gaussian_shaped_labels(sigma, sz)
%GAUSSIAN_SHAPED_LABELS
%   Gaussian-shaped labels for all shifts of a sample.
%
%   LABELS = GAUSSIAN_SHAPED_LABELS(SIGMA, SZ)
%   Creates an array of labels (regression targets) for all shifts of a
%   sample of dimensions SZ. The output will have size SZ, representing
%   one label for each possible shift. The labels will be Gaussian-shaped,
%   with the peak at 0-shift (top-left element of the array), decaying
%   as the distance increases, and wrapping around at the borders.
%   The Gaussian function has spatial bandwidth SIGMA.
%为SZ维度样本的所有位移创建一个标签数组（回归目标）。 输出的尺寸为SZ，代表每个可能的移位的一个标签。 
%标签将是高斯形状的，其峰值为0移位（阵列的左上角元素），随距离增加而衰减，并在边界环绕。高斯函数具有空间带宽SIGMA。
%   Joao F. Henriques, 2014
%   http://www.isr.uc.pt/~henriques/


% 	%as a simple example, the limit sigma = 0 would be a Dirac delta,
% 	%instead of a Gaussian:
% 	labels = zeros(sz(1:2));  %labels for all shifted samples
% 	labels(1,1) = magnitude;  %label for 0-shift (original sample)


%evaluate a Gaussian with the peak at the center element用中心元素处的峰值
[rs, cs] = ndgrid((1:sz(1)) - floor(sz(1)/2), (1:sz(2)) - floor(sz(2)/2));
labels = exp(-0.5 / sigma^2 * (rs.^2 + cs.^2));

%move the peak to the top-left, with wrap-around将峰顶移至左上角，并进行环绕
labels = circshift(labels, -floor(sz(1:2) / 2) + 1);
% B= circshift(A,K,m);
% 输入参数：A表示待移位的矢量或矩阵；
% K表示所移位数，可以是数字，也可以是二维数组，若是数字则可以和m协同作用来决定是行移位还是列移
% m当K是数字时，m用来决定是行移位还是列移位。默认m是1，当m=1时表示列移位，当m=2时表示行移       
% 当K=1时，分别对每列元素向前移1位   -1   向后   直接令m=2则表示行移位  K= [col，row]，其中col表示列位移，row表示行位移

%sanity check: make sure it's really at top-left
assert(labels(1,1) == 1)

end

