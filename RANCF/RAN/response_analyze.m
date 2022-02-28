
function [ state ,score ] = response_analyze( response )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
global ran
global enableGPU

if isempty(ran)
    initial_RAN();
end

response = single(response);
response = imResample(response, ran.meta.inputSize(1:2));

ran.layers{end}.type = 'softmax';

if enableGPU
    response = gpuArray(response); 
    ran = vl_simplenn_move(ran,'gpu');
end

res = vl_simplenn(ran,response);
scores = squeeze(gather(res(end).x));
[~, state] = max(scores);
score = scores(2);

end

