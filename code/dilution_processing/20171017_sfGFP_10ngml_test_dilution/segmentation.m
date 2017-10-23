% Segment the autofluorescence and the original images using SuperSegger
% in MATLAB.

% Define the experiment parameters.
DATE = '20171017';
BASENAME = 'sfGFP_10ngmL';
samples = {'autofluorescence', 'growth'};
CONST = loadConstants('60XCaulob');
CONST.trackFoci.numSpots = 0;
CONST.align.ALIGN_FLAG = 1;
CONST.trackOpti.REMOVE_STRAY = 1;
cleanFlag = 1;
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