% GET_FEATURES: Extracting hierachical convolutional features
%indLayers = [37, 28, 19];
function feat = get_features(im, cos_window, layers)
%feat  = get_features(patch, cos_window, indLayers);%indLayers = [37, 28, 19];
global net
global enableGPU

if isempty(net)
    initial_net();
end

sz_window = size(cos_window);

% Preprocessing
img = single(im);        % note: [0, 255] range,和uint8范围一样  single占4个字节 比double类型省一半内存
img = imResample(img, net.meta.normalization.imageSize(1:2));%调整图片大小



average=net.meta.normalization.averageImage;

if numel(average)==3
    average=reshape(average,1,1,3);
end

img = bsxfun(@minus, img, average);

if enableGPU, img = gpuArray(img); net = vl_simplenn_move(net,'gpu');end

% Run the CNN

res = vl_simplenn(net,img);

% Initialize feature maps
feat = cell(length(layers), 1);%length是求某一矩阵所有维的最大长度，初始化feat=cell(3,1);

for ii = 1:length(layers)
    
    % Resize to sz_window
    if enableGPU
        x = gather(res(layers(ii)).x); 
    else
        x = res(layers(ii)).x;
    end
    
    x = imResample(x, sz_window(1:2));
    
    % windowing technique
    if ~isempty(cos_window)
        x = bsxfun(@times, x, cos_window);
    end
    
    feat{ii}=x;%feat返回卷积特征层，共3层
end

end
