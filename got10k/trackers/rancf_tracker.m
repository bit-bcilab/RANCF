
function [boxes, times] = rancf_tracker(img_files, box, visualize)
    frames = length(img_files);
    target_sz = [box(1,4), box(1,3)];
    pos = [box(1,2), box(1,1)] + floor(target_sz/2);
    
    % 参数预设置
    padding = struct('generic', 1.8, 'large', 1, 'height', 0.4);
    hog_padding = struct('generic', 1.5, 'large', 1.5, 'height', 1.3);

    lambda = 1e-4;              % Regularization parameter (see Eqn 3 in our paper)正则化参数（见本文中的公式3）
    output_sigma_factor = 0.1;  % Spatial bandwidth (proportional to the target size)空间带宽（与目标大小成比例）

    interp_factor = 0.012;       % Model learning rate (see Eqn 6a, 6b)模型学习率（见公式6a，6b）
    cell_size = 4;              % Spatial cell size空间单元格大小
    
    % 网络选择、卷积层选择、权重设置
    global enableGPU;
    enableGPU = true;

    [boxes, times] = CoMCCF_tracker(img_files, pos, target_sz, ...
        padding, lambda, output_sigma_factor, interp_factor, cell_size, hog_padding);
end


function [rects, time] = CoMCCF_tracker(frames, pos, target_sz, ...
    padding, lambda, output_sigma_factor, interp_factor, cell_size, hog_padding)


% ================================================================================
% Environment setting
% ================================================================================

indLayers = [37, 28, 19];   % The CNN layers Conv5-4, Conv4-4, and Conv3-4 in VGG Net
nweights  = [1, 0.5, 0.25]; % Weights for combining correlation filter responses 融合权重
numLayers = length(indLayers);%3层

% 图片的尺寸
im = imread([frames{1}]);
if ismatrix(im)
    im = cat(3, im, im, im);
end
im_sz = size(im);
output_sigma = sqrt(prod(target_sz)) * output_sigma_factor / cell_size;%3.0873

window_sz = get_search_window(target_sz, im_sz, padding);
l1_patch_num = floor(window_sz/ cell_size);%floor 高斯取整  不大于的整数  [62 61]
yf = fft2(gaussian_shaped_labels(output_sigma, l1_patch_num)); % 62x61double
cos_window = hann(size(yf,1)) * hann(size(yf,2))';

hog_window_sz = get_search_window_MCCF(target_sz, im_sz, hog_padding);
hog_l1_patch_num = floor(hog_window_sz/ cell_size);%floor 高斯取整  不大于的整数  [62 61]
hog_yf = fft2(gaussian_shaped_labels(output_sigma, hog_l1_patch_num)); % 62x61double
hog_cos_window = hann(size(hog_yf,1)) * hann(size(hog_yf,2))';
                

% Initialize variables for calculating FPS and distance precision
time = 0;
rects = zeros(numel(frames), 4);
nweights  = reshape(nweights,1,1,[]);


% ================================================================================
% Initialization
% ================================================================================

model_xf = cell(1, numLayers);
model_af = cell(1, numLayers);
score = zeros(numel(frames), 1);

hog_score = zeros(numel(frames), 1);
dist = zeros(numel(frames), 1);
diag = sqrt(target_sz(1)^2+target_sz(2)^2);
hog_power = 1;
help_num = 0;

is_occ = 0;

feat  = extractFeature(im, pos, window_sz, cos_window, indLayers);
[model_xf, model_af] = updateModel(feat, yf, interp_factor, lambda,...
                                 'init', model_xf, model_af);

hog_xf = cell(1, 1);
hog_af = cell(1, 1);
hog_feat = hog_extractFeature(im, pos, hog_window_sz, hog_cos_window, cell_size);
[hog_xf, hog_af] = updateModel(hog_feat, hog_yf, interp_factor, lambda,...
                                 'init', hog_xf, hog_af);
                             
hog_pos = pos;
                             
current_scale_factor=1;
init_scale_para(im, target_sz, pos);                             
target_sz_t=target_sz*current_scale_factor;
rects(1,:) = [pos([2,1]) - target_sz_t([2,1])/2, target_sz_t([2,1])];

% ================================================================================
% Start tracking
% ================================================================================

for frame = 2:numel(frames)
    
    im = imread([ frames{frame}]); % Load the image at the current frame
    if ismatrix(im)
        im = cat(3, im, im, im);
    end
    is_hog = 0;
    tic();
    
    % ================================================================================
    % Predicting the object position from the learned object model
    % ================================================================================
        
    % Extracting hierarchical convolutional features提取分层卷积特征
    feat = extractFeature(im, pos, window_sz, cos_window, indLayers);
    [pos_, response]  = predictPosition(feat, pos, indLayers, nweights, cell_size,...
                                                 l1_patch_num, model_xf, model_af);                                        
    [~ ,score(frame)] = response_analyze(response);
    
    if frame == 50
       hog_power = mean(hog_score(1:50)); 
    end
    
    if hog_power > 0.45 && help_num < 1
        hog_feat = hog_extractFeature(im, hog_pos, hog_window_sz, hog_cos_window, cell_size);
        [hog_pos_, hog_response]  = predictPosition(hog_feat, hog_pos, [0], [1], cell_size,...
            hog_l1_patch_num, hog_xf, hog_af);
        hog_response = hog_response * 1.65;
        [~ ,hog_score(frame)] = response_analyze(hog_response);
        
        dist(frame) = sqrt((hog_pos_(1) - pos_(1))^2+(hog_pos_(2) - pos_(2))^2);
    end
    
    % 遮挡判断
    if (is_occ == 0) && (frame > 10) && (score(frame) < 0.3) && ((score(frame-3)-score(frame)) > 0.38) && (var(score(frame-10:frame-3)) < 0.006)
        is_occ = 1;
        occ_pos = pos;
        occ_xf = model_xf;
        occ_af = model_af;    
    end
    
    % 结果中记录位置以主分类器位置为主
    pos = pos_;
    
    % 遮挡时，分别在当前位置和消失位置上进行定位，哪个位置找到就在哪个位置继续跟踪，遮挡时模型停止更新
    if is_occ        

        % 遮挡分类器永远在目标消失的位置附近进行检测
        occ_feat = extractFeature(im, occ_pos, window_sz, cos_window, indLayers);
        
        [~, occ_response]  = predictPosition(occ_feat, occ_pos, indLayers, nweights, cell_size,...
                                                      l1_patch_num, occ_xf, occ_af);
        [~ ,occ_score] = response_analyze(occ_response);
        
        % 主分类器响应大于0.3时，遮挡停止，继续正常跟踪流程
        if score(frame) > 0.3 
            is_occ = 0;
            hog_pos = occ_pos;
            % 遮挡分类器响应大于0.3时，遮挡停止，从消失位置继续跟踪
        elseif occ_score > 0.3
            is_occ = 0;
            pos = occ_pos;
            hog_pos = occ_pos;
        end
   
    else
        dist_hog = sqrt((hog_pos_(1) - hog_pos(1))^2 + (hog_pos_(2) - hog_pos(2))^2);
        if frame > 10
            if dist(frame)> diag/2 && hog_power > 0.4 && help_num < 1
                if mean(score(frame-1:frame)) >  mean(hog_score(frame-1:frame)) && score(frame) > hog_score(frame)
                    hog_pos_ = pos_;
                    help_num = help_num + 1;
                else
                    if hog_score(frame) > 0.3 && dist_hog < diag/2 && var(hog_score(frame-10:frame)) < 0.01
                        pos = hog_pos_;
                    end
                end
            end
        end

        if hog_power > 0.45 && help_num < 1
            hog_pos = hog_pos_;
            hog_feat = hog_extractFeature(im, hog_pos, hog_window_sz, hog_cos_window, cell_size);
            hog_if = 0.01;
            [hog_xf, hog_af] = updateModel(hog_feat, hog_yf, 2*hog_if, lambda,...
                'main', hog_xf, hog_af);
        end
            
        if (score(frame) > 0.2) || is_hog
            
            if frame > 10 && sum(score(frame-4:frame)>0.85) == 5
                if sum(score(frame-7:frame-5)>0.85) == 3
                    if mod(frame,3) == 0
                        upd = 1;
                        addit = 2.5;
                    else
                        upd = 0;
                        addit = 1;
                    end
                else
                    if mod(frame,2) == 0
                        upd = 1;
                        addit = 2;
                    else
                        upd = 0;
                        addit = 1;
                    end
                end
            else
                upd = 1;
                addit = 1;
            end     

            if upd
                
                feat = extractFeature(im, pos, window_sz, cos_window, indLayers);
                
                if frame < 250
                    interp_factor = (1/100) * (1 / (1 + exp(-score(frame))) + 0.325) * addit;
                else
                    interp_factor = (1/100) * (1 / (1 + exp(-score(frame)) + 0.35))  * addit;
                end
                
                
                [model_xf, model_af] = updateModel(feat, yf, interp_factor, lambda,...
                    'main', model_xf, model_af);
                
            
            end
        end
        
    end
       
    % Scale estimation 比例估算
    current_scale_factor = estimate_scale( im, pos, current_scale_factor);   
    
    % ================================================================================
    % Save predicted position and timing
    % ================================================================================

    target_sz_t=target_sz*current_scale_factor;
    rects(frame,:) = [pos([2,1]) - target_sz_t([2,1])/2, target_sz_t([2,1])];

    time = time + toc();
    
end

end


function [pos,response] = predictPosition(feat, pos, indLayers, nweights, cell_size, l1_patch_num, ...
    model_xf, model_alphaf)
    %cell(1,3)    %cell(1,3)
% ================================================================================
% Compute correlation filter responses at each layer
% ================================================================================
res_layer = zeros([l1_patch_num, length(indLayers)]);

for ii = 1 : length(indLayers)%第5层，到4层，到3层
    zf = fft2(feat{ii});%对卷积特征层傅里叶变换
    kzf=sum(zf .* conj(model_xf{ii}), 3) / numel(zf);
    %y=conj(x)函数计算复数x的共轭值。输出结果y的维数跟输入x的维数一致
    %n=numel(A)该语句返回数组中元素的总数
    %A为三通道图像，则sum(A,3)运算后的值为每个通道对应位置的值各自相加，比如在位置p三通道像素值分别为r,g,b,则在p位置运算后的值为r+g+b
    temp= real(fftshift(ifft2(model_alphaf{ii} .* kzf)));  %equation for fast detection快速检测方程
    %fftshift将零频率的分量移到频谱的中心语法  
    res_layer(:,:,ii)=temp/max(temp(:));%temp(1,2,3)就是公式中的f(5,4,3)层,res这里有归一化操作
end
% Combine responses from multiple layers (see Eqn. 5)组合来自多个层的响应（参见公式5）
response = sum(bsxfun(@times, res_layer, nweights), 3);


%C=bsxfun(fun,A,B)：两个数组间元素逐个计算  @times   Array multiply
%nweights  = [1, 0.5, 0.25]; 对应5 4 3层

% ================================================================================
% Find target location
% ================================================================================
% Target location is at the maximum response. we must take into
% account the fact that, if the target doesn't move, the peak
% will appear at the top-left corner, not at the center (this is
% discussed in the KCF paper). The responses wrap around cyclically.
[vert_delta, horiz_delta] = find(response == max(response(:)), 1);
vert_delta  = vert_delta  - floor(size(zf,1)/2);
horiz_delta = horiz_delta - floor(size(zf,2)/2);

% Map the position to the image space
pos = pos + cell_size * [vert_delta - 1, horiz_delta - 1];


end


function [model_xf, model_alphaf] = updateModel(feat, yf, interp_factor, lambda, state, ...
    model_xf, model_alphaf)

numLayers = length(feat);

% ================================================================================
% Initialization
% ================================================================================
xf       = cell(1, numLayers);
alphaf   = cell(1, numLayers);

% ================================================================================
% Model update
% ================================================================================
for ii=1 : numLayers % 5 4 3 层
    xf{ii} = fft2(feat{ii});%对feat卷积特征FFT
    kf = sum(xf{ii} .* conj(xf{ii}), 3) / numel(xf{ii});
    alphaf{ii} = yf./ (kf+ lambda);   % Fast training
end

switch state
    % 第一帧的初始化
    case 'init'
        for ii=1:numLayers
            model_alphaf{ii} = alphaf{ii};
            model_xf{ii}     = xf{ii};
        end
        % 主分类器的更新
    case 'main'
        for ii=1:numLayers
            model_alphaf{ii} = (1 - interp_factor) * model_alphaf{ii} + interp_factor * alphaf{ii};
            model_xf{ii}     = (1 - interp_factor) * model_xf{ii}     + interp_factor * xf{ii};
        end
end

% Online model update using learning rate interp_factor




end

function feat  = extractFeature(im, pos, window_sz, cos_window, indLayers)

% Get the search window from previous detection从以前的检测中获取搜索窗口
patch = get_subwindow(im, pos, window_sz);%取了输入图像的目标框内的矩阵  取了window_sz，以pos为中心点
% Extracting hierarchical convolutional features   patch  250x244x3
feat  = get_features(patch, cos_window, indLayers);%indLayers = [37, 28, 19];

end

function feat  = hog_extractFeature(im, pos, window_sz, cos_window, cell_size)

% Get the search window from previous detection从以前的检测中获取搜索窗口
patch = get_subwindow(im, pos, window_sz);%取了输入图像的目标框内的矩阵  取了window_sz，以pos为中心点


feat = cell(1,1);
feat{1} = double(fhog(single(patch) / 255, cell_size, 9));
feat{1}(:,:,end) = [];  %remove all-zeros channel ("truncation feature")
feat{1} = bsxfun(@times, feat{1}, cos_window);


end

