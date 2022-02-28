setup_got10k;

% tracker_fn = @identity_tracker;
% experiment = ExperimentGOT10k('F:\DataBase\GOT-10k', 'test');
% experiment.run('IdentityTracker', tracker_fn, false);
% experiment.report({'IdentityTracker'});

tracker_fn = @rancf_tracker;
experiment = ExperimentLaSOT('K:\LaSOT\LaSOTBenchmark');
experiment.run('RANCF', tracker_fn, false);

% tracker_fn = @rancf_tracker;
% experiment = ExperimentOTB('F:\DataBase\Benchmark');
% experiment.run('RANCF', tracker_fn, false);
% 
% tracker_fn = @rancf_tracker;
% experiment = ExperimentGOT10k('F:\DataBase\GOT-10k', 'test');
% experiment.run('RANCF', tracker_fn, false);
% experiment.report({'RANCF'});
