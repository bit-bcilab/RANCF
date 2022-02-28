function [img_files, pos, target_sz, ground_truth, video_path, rect] = load_video_info(base_path, video)
%LOAD_VIDEO_INFO
%   Loads all the relevant information for the video in the given path:
%   the list of image files (cell array of strings), initial position图像文件列表（字符串的单元格数组），初始位置
%   (1x2), target size (1x2), the ground truth information for precision（1x2），目标大小（1x2），精确度的地面实况信息
%   calculations (Nx2, for N frames), and the path where the images are计算（Nx2，N帧）以及图像所在的路径
%   located. The ordering of coordinates and sizes is always [y, x].位于。 坐标和大小的排序总是[y，x]。
%
%   Joao F. Henriques, 2014
%   http://www.isr.uc.pt/~henriques/


	%see if there's a suffix, specifying one of multiple targets, for
	%example the dot and number in 'Jogging.1' or 'Jogging.2'.
	if numel(video) >= 2 && video(end-1) == '.' && ~isnan(str2double(video(end)))
		suffix = video(end-1:end);  %remember the suffix
		video = video(1:end-2);  %remove it from the video name
	else
		suffix = '';
	end

	%full path to the video's files
	if base_path(end) ~= '/' && base_path(end) ~= '\'
		base_path(end+1) = '/';
	end
	video_path = [base_path video '/'];

	%try to load ground truth from text file (Benchmark's format)尝试从文本文件加载基本事实（基准格式）
% 	filename = [video_path 'groundtruth' suffix '.txt'];
    filename = [video_path 'groundtruth_rect' suffix '.txt'];
    f = fopen(filename);
    
    if f == -1
        filename = [video_path 'groundtruth' suffix '.txt'];
        f = fopen(filename);
    end
    
    % 2018.8.21 22:00 修改了load_video_info函数，使其可以顺利读入未标注的序列(1)
    sqe_signed = 1;%判断序列是否被标注过
    
	%the format is [x, y, width, height]格式
    if f == -1
        sqe_signed = 0;
        target_sz = [0,0];
        pos = [0,0];
    else
        try
            ground_truth = load(filename);
        catch  %try different format (no commas)
            frewind(f);
            ground_truth = textscan(f, '%f %f %f %f ');
        end
        fclose(f);
        
        [~,D] = size(ground_truth);
        if D == 8
            ground_truth(:,1:2) = ground_truth(:,3:4);
            ground_truth(:,3:4) = ground_truth(:,7:8) - ground_truth(:,3:4);
            ground_truth(:,5:8) = [];
        end
        
        %set initial position and size
        rect = ground_truth;
        target_sz = [ground_truth(1,4), ground_truth(1,3)];%目标大小(高，宽）
        pos = [ground_truth(1,2), ground_truth(1,1)] + floor(target_sz/2);%中心点坐标（y,x) (行，列）
        
        if size(ground_truth,1) == 1 %size(x,dim) dim==1,返回行数
            %we have ground truth for the first frame only (initial position)
            ground_truth = [];
        else
            %store positions instead of boxes存储位置而不是框
            ground_truth = ground_truth(:,[2,1]) + ground_truth(:,[4,3]) / 2;
        end

    end
	
	
	%for these sequences, we must limit ourselves to a range of frames.对于这些序列，我们必须将自己限制在一系列帧中。
	%for all others, we just load all png/jpg files in the folder.所有其他的，我们只加载文件夹中的所有png / jpg文件。
	frames = {'David', 300, 770;
			  'Football1', 1, 74;
			  'Freeman3', 1, 460;
			  'Freeman4', 1, 283;
              'Diving', 1, 215};
	idx = find(strcmpi(video, frames(:,1)));
	
	if isempty(idx)
		%general case, just list all images
		img_files = dir([video_path '*.png']);
        if isempty(img_files)
            img_files = dir([video_path 'img/' '*.png']);
            if isempty(img_files)
                img_files = dir([video_path '*.jpg']);
                if isempty(img_files)
                    img_files = dir([video_path 'img/' '*.jpg']);
                    assert(~isempty(img_files), 'No image files to load.');
                    video_path = [video_path 'img/'];
                end
            else
                video_path = [video_path 'img/'];
            end
        end
		img_files = sort({img_files.name});
        
        % 未标注序列的读取
        if sqe_signed == 0
            ground_truth = zeros(numel(img_files),2);
            rect = zeros(numel(img_files),4);
            [pos, target_sz ] = read_notsigned( video_path,img_files);
        end
	else
		%list specified frames. try png first, then jpg.
        video_path = [video_path 'img/'];
		if exist(sprintf('%s%04i.png', video_path, frames{idx,2}), 'file')
			img_files = num2str((frames{idx,2} : frames{idx,3})', '%04i.png');
			
		elseif exist(sprintf('%s%04i.jpg', video_path, frames{idx,2}), 'file')
			img_files = num2str((frames{idx,2} : frames{idx,3})', '%04i.jpg');
			
		else
			error('No image files to load.')
		end
		
		img_files = cellstr(img_files);
        
	end
	
end

