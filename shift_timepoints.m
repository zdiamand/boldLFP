function shifted_timepoints = shift_timepoints(timepoints,offset)
    shifted_timepoints = timepoints;
    shifted_timepoints(timepoints ~= 0) = timepoints(timepoints ~= 0) - timepoints(1) + offset;
end