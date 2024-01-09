function [conflicting_bids_curr] = conflicting_bids_with_curr_alloc_modified(curr_alloc,conflicting_bids_array,PARTICLE_WIDTH_x,PARTICLE_WIDTH_y,prev_conflicting_bids)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% - INPUT
%%%%%%          curr_alloc - an array of bids in curr allocation in the
%%%%%%          first column, among with their requested goods in columns
%%%%%%          3:end.
%%%%%%          conflicting_bids_array - each row is a good num, in cols
%%%%%%          2:end there're the bids requesting that good.
%%%%%%          for example = all the bids wanting good number 7, are located in
%%%%%%          conflicting_bids_array(7,2:end).
%%%%%%          prev_conflicting_bids - in order to avoid the calculation
%%%%%%          of the same thing all of the time, we get as input the bids
%%%%%%          that conflict with the bids in curr_alloc(1:end-1,:) -
%%%%%%          because we already computed those bids in the prevoius run
%%%%%%          and we only need to update the new conflicted bids due to
%%%%%%          the bid we've added in curr_alloc(end)
%%%%%% - OUTPUT
%%%%%%          conflicting_bids - array of bid nums conflicting with that
%%%%%%          curr_allocation - including the curr allocation itself.
%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%always run on one bid

goods_requested_by_curr_bid = curr_alloc(end, 3:end);
% goods_requested_by_curr_bid = goods_requested_by_curr_bid(isfinite(goods_requested_by_curr_bid));

conflicting_bids_indices_2 = conflicting_bids_array(goods_requested_by_curr_bid(isfinite(goods_requested_by_curr_bid)), 2:end);
conflicting_bids_indices = conflicting_bids_indices_2(:);

% Merge conflicting bids
conflicting_bids_curr = (unique([prev_conflicting_bids; conflicting_bids_indices]));

end