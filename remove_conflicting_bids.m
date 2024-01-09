function updated_bids_array = remove_conflicting_bids(bids_array,conflicting_bids_with_curr_alloc)
updated_bids_array = bids_array;
matching_rows = ismember(bids_array(:, 1), conflicting_bids_with_curr_alloc);
updated_bids_array(matching_rows, :) = NaN;
updated_bids_array = updated_bids_array(~matching_rows, :);
end