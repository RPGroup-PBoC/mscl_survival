% Segment the autofluorescence and the original images using SuperSegger
% in MATLAB

% Define the experiment parameters.
DATE = '20171129';
BASENAME = 'sfGFP_10ngmL';
samples = {'auto', 'growth'};
CONST = loadConstants('60XCaulob');
CONST.parallel.PARALLEL_FLAG = 0;
CONST.trackFoci.numSpots = 0;
CONST.align.ALIGN_FLAG = 1;
CONST.trackOpti.REMOVE_STRAY = 1;
cleanFlag = 0;
for i=1:length(samples)
    statement = ['Beginning segmentaton ', num2str(i), ' out of ',...
        num2str(length(samples))];
    disp(statement)
    % Define the data directory.
    directory = ['../../../data/images/', DATE, '_', BASENAME,...
        '_dilution/', samples{i}, '/'];

    % Begin the segmentation.
    BatchSuperSeggerOpti(directory, 1, cleanFlag, CONST);
end

disp('Finished!');
