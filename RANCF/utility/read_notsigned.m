function [ pos, target_sz ] = read_notsigned( video_path,img_files )
%read_notsigned
%   未标注序列的读取

    im = imread([video_path img_files{1}]);
    figure(1)
    imshow(im);
    [temp, rect] = imcrop(im);
    [target_sz(1), target_sz(2),~] = size(temp);
    pos = round([rect(2), rect(1)] + target_sz/2);
    close(figure(1));
    
end

