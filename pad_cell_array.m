function padded_cell_array = pad_cell_array(cell_array)
    % Define the max length of the double arrays
    max_len = max(cellfun(@length, cell_array));

    % Preallocate a new cell array with the padded arrays
    padded_cell_array = cell(1, 300);

    % Loop through each entry in the original cell array
    for i = 1:300
        d = cell_array{i};
        % Check if d is a 1 by d or d by 1 array
        if size(d, 1) == 1
            % Transpose the array to make it a d x 1 double
            d = d.';
        end
        % Pad the current entry with zeros so it has the max length
        padded_d = [d; zeros(max_len-length(d), 1)];
        % Transpose the padded array to make it a 1xd double
        padded_d = padded_d.';
        % Store the padded array in the new cell array
        padded_cell_array{i} = padded_d;
    end
end
