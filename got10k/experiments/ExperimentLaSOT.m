classdef ExperimentLaSOT
    
    properties
        dataset
        result_dir
        report_dir
        nbins_iou
        repetitions
        nbins_ce
    end

    methods

        function obj = ExperimentLaSOT(root_dir, version, result_dir, report_dir)
            if nargin < 2
                version = 2015;
            end
            if nargin < 3
                result_dir = 'results';
            end
            if nargin < 4
                report_dir = 'reports';
            end

            obj.dataset = LaSOT(root_dir, version, false);
            obj.result_dir = fullfile(result_dir, 'LaSOT');
            obj.report_dir = fullfile(report_dir, 'LaSOT');
            % as nbins_iou increases, the success score
            % converges to average overlap (AO)
            obj.nbins_iou = 101;
            obj.nbins_ce = 51;
        end

        function obj = run(obj, tracker_name, tracker_fn, visualize)
            if nargin < 3
                visualize = false;
            end
            fprintf('Running tracker %s on OTB...\n', tracker_name);
            
            record_dir = fullfile(obj.result_dir, sprintf('%s_tracking_result', tracker_name));
            mkdir(record_dir)
            % loop over the complete dataset
            for s = 1:length(obj.dataset)
                seq_name = obj.dataset.seq_names{s};
                fprintf('--Sequence %d/%d: %s\n', s, length(obj.dataset), seq_name);

                % tracking loop
                [img_files, anno] = obj.dataset(s);
                [rects, speed_fps] = tracker_fn(img_files, anno(1, :), visualize);
                assert(size(rects, 1) == size(anno, 1));

                % record results
                record_file = fullfile(record_dir, sprintf('%s.txt', seq_name));
                writematrix(rects, record_file)
                dlmwrite(record_file, rects, 'delimiter', ',', 'precision', '%.3f');
                
%                 obj.record(tracker_fn, seq_name, rects, speed_fps);
            end
        end

        function obj = record(obj, tracker_name, seq_name, rects, speed_fps)
        end

    end

end
