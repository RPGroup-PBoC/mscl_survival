function status = dilutionSegmentation(dirName, const_file)
% DILUTIONSEGMENTATION This function segments E. coli cells growing on an M9 +
% 0.5% glucose agarose pad and extracts fluorescence information.

% Load the constants file and set parameters.
CONST = loadConstants(const_file);
CONST.align.ALIGN_FLAG = 1;

% Launch the segmentation on the provided directory.
BatchSuperSeggerOpti(dirName, 1, 1, CONST);

status = 'Finished';
end
