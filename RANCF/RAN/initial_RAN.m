
function initial_RAN()
% INITIAL_NET: Loading VGG-Net-19

global ran;
ran = load('RAN.mat');
ran =ran.net;

% Switch to GPU mode
global enableGPU;
if enableGPU
    ran = vl_simplenn_move(ran, 'gpu');
end

ran=vl_simplenn_tidy(ran);

end