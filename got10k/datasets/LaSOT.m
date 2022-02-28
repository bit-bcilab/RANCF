classdef LaSOT

    properties
        otb13_seqs
        otb15_seqs
        tb50_seqs
        tb100_seqs
        root_dir
        version
        seq_names
        seq_dirs
        anno_files
    end
    
    methods

        function obj = LaSOT(root_dir, version, download)
            if nargin < 2
                % version has to be one of 2013, 2015,
                % 'otb2013', 'otb2015', 'tb50' and 'tb100'
                version = 2015;
            end
            if nargin < 3
                download = true;
            end
            
            obj.root_dir = root_dir;
            obj.version = version;
            if download
                obj.download(root_dir, version);
            end
            obj.check_integrity(root_dir, version);

            valid_seqs = obj.get_seqs(version);
            obj.anno_files = {};
            for s = 1:length(valid_seqs)
                valid_seq = valid_seqs{s};
                father_dir = valid_seq(1:strfind(valid_seq, '-')-1);
                files = dir(fullfile(root_dir, father_dir, valid_seq, 'groundtruth*.txt'));
                obj.anno_files = [obj.anno_files fullfile({files.folder}, {files.name})];
            end
            % remove empty annotation files
            % (e.g., groundtruth_rect.1.txt of Human4)
            obj.anno_files = obj.filter_files(obj.anno_files);
            obj.seq_dirs = {};
            obj.seq_names = {};
            for f = 1:length(obj.anno_files)
                obj.seq_dirs{f} = fileparts(obj.anno_files{f});
                [~, obj.seq_names{f}, ~] = fileparts(obj.seq_dirs{f});
            end
            % rename repeated sequence names
            % (e.g., Jogging and Skating2)
            obj.seq_names = obj.rename_seqs(obj.seq_names);
        end

        function varargout = subsref(obj, s)
            switch s(1).type
                case '.'
                    varargout{1} = builtin('subsref', obj, s);
                case {'()', '{}'}
                    s.type = '{}';
                    seq_dir = builtin('subsref', obj.seq_dirs, s);
                    img_files = dir(fullfile(seq_dir, 'img/*.jpg'));
                    img_files = sort(fullfile(...
                        {img_files.folder}, {img_files.name}));

                    % special sequences
                    % (visit http://cvlab.hanyang.ac.kr/tracker_benchmark/index.html for detail)
                    [~, seq_name, ~] = fileparts(seq_dir);
                    if strcmpi(seq_name, 'david')
                        img_files = img_files(300:770);
                    elseif strcmpi(seq_name, 'football1')
                        img_files = img_files(1:74);
                    elseif strcmpi(seq_name, 'freeman3')
                        img_files = img_files(1:460);
                    elseif strcmpi(seq_name, 'freeman4')
                        img_files = img_files(1:283);
                    elseif strcmpi(seq_name, 'diving')
                        img_files = img_files(1:215);
                    end

                    anno_file = builtin('subsref', obj.anno_files, s);
                    anno = dlmread(anno_file);
                    assert(length(img_files) == size(anno, 1));
                    varargout{1} = img_files;
                    varargout{2} = anno;
            end
        end

        function seq_num = length(obj)
            seq_num = length(obj.seq_names);
        end

        function seq_names = get_seqs(obj, version)
            seq_names = config_sequence('test_set');
        end

        function filtered_files = filter_files(obj, filenames)
            filtered_files = {};
            for f = 1:length(filenames)
                content = fileread(filenames{f});
                content = strip(content);
                if strcmp(content, '')
                    fprintf('warning: %s is empty\n', filenames{f});
                else
                    filtered_files = [filtered_files filenames{f}];
                end
            end
        end

        function renamed_seqs = rename_seqs(obj, seq_names)
            % in case some sequences may have multiple targets
            renamed_seqs = {};
            for s = 1:length(seq_names)
                if sum(strcmp(seq_names{s}, seq_names)) == 1
                    renamed_seqs = [renamed_seqs seq_names{s}];
                else
                    ind = sum(strcmp(seq_names{s}, seq_names(1:s)));
                    renamed_seqs = [renamed_seqs sprintf('%s.%d', seq_names{s}, ind)];
                end
            end
        end

        function root_dir = download(obj, root_dir, version)
            seq_names = obj.get_seqs(version);

            if ~exist(root_dir, 'dir')
                mkdir(root_dir);
            else
                downloaded = true;
                for s = 1:length(seq_names)
                    if ~exist(fullfile(root_dir, seq_names{s}), 'dir')
                        downloaded = false;
                        break;
                    end
                end
                if downloaded
                    disp('Files already downloaded.');
                    return;
                end
            end

            url_fmt = 'http://cvlab.hanyang.ac.kr/tracker_benchmark/seq/%s.zip';
            for s = 1:length(seq_names)
                seq_dir = fullfile(root_dir, seq_names{s});
                if exist(seq_dir, 'dir')
                    continue;
                end
                url = sprintf(url_fmt, seq_names{s});
                zip_file = fullfile(root_dir, [seq_names{s} '.zip']);
                fprintf('Downloading to %s...\n', zip_file);
                urlwrite(url, zip_file);
                fprintf('Extracting to %s...\n', root_dir);
                unzip(zip_file, root_dir);
            end
        end

        function check_integrity(obj, root_dir, version)
            seq_names = obj.get_seqs(version);

            if exist(root_dir, 'dir')
                % check each sequence folder
                for s = 1:length(seq_names)
                    seq_dir = fullfile(root_dir, seq_names{s});
                    if ~exist(seq_dir, 'dir')
                        fprintf('Warning: sequence %s not exist.', seq_names{s});
                    end
                end
            else
                % dataset not exist
                error(['Dataset not found or corrupted. '...
                       'You can set download to true to download it.']);
            end
        end

    end

end
