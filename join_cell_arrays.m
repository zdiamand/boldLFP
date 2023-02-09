function joined_cell_array = join_cell_arrays(cell_array1, cell_array2)
    % Preallocate a new cell array with the joined arrays
    joined_cell_array = cell(1, 300);

    % Loop through each entry in the two original cell arrays
    for i = 1:300
        d1 = cell_array1{i};
        d2 = cell_array2{i};
        % Join the current entries from the two cells
        joined_d = [d1; d2];
        % Store the joined array in the new cell array
        joined_cell_array{i} = joined_d;
    end
end
