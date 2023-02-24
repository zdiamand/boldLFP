function joined_cell_array = join_cell_arrays(varargin)
    % Get the number of cell arrays
    num_arrays = nargin;
    % Find the length of the cell arrays
    array_length = 300;
    
    % Preallocate a new cell array with the joined arrays
    joined_cell_array = cell(1, array_length);
    
    % Loop through each entry in the original cell arrays
    for i = 1:array_length
        joined_d = [];
        for j = 1:num_arrays
            current_array = varargin{j};
            d = current_array{i};
            % Join the current entries from the cell arrays
            joined_d = [joined_d; d];
        end
        % Store the joined array in the new cell array
        joined_cell_array{i} = joined_d;
    end
end