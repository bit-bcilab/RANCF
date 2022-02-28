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
img = single(im);        % note: [0, 255] range,��uint8��Χһ��  singleռ4���ֽ� ��double����ʡһ���ڴ�
img = imResample(img, net.meta.normalization.imageSize(1:2));%����ͼƬ��С



average=net.meta.normalization.averageImage;

if numel(average)==3
    average=reshape(average,1,1,3);
end

img = bsxfun(@minus, img, average);

if enableGPU, img = gpuArray(img); net = vl_simplenn_move(net,'gpu');end

% Run the CNN

res = vl_simplenn(net,img);

% Initialize feature maps
feat = cell(length(layers), 1);%length����ĳһ��������ά����󳤶ȣ���ʼ��feat=cell(3,1);

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
    
    feat{ii}=x;%feat���ؾ�������㣬��3��
end

end
