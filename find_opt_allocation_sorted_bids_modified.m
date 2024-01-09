function [opt_comb_allocation, opt_revenue,run_time,naive_allocation, naive_revenue, full_allocations_explored_with_th, num_of_pruning_conditions_calc] = find_opt_allocation_sorted_bids_modified(bids_array,NUM_OF_PARTICLES,PARTICLE_WIDTH_x,PARTICLE_WIDTH_y)
tic;
%% sort bids by max vals

sorted_bids_array = sortrows(bids_array,2,"descend");
num_of_bids = length(bids_array(:,1));
conflicting_bids_array = conflicting_bids(bids_array,PARTICLE_WIDTH_x,PARTICLE_WIDTH_y);
% sorted_bids_array = sortrows(bids_array,2,"descend");
opt_revenue = 0;
%%
max_vals = maxk(sorted_bids_array(:,2),NUM_OF_PARTICLES);

bids_array_no_i         =   cell(NUM_OF_PARTICLES,1);
bids_left_num           =   cell(NUM_OF_PARTICLES,1);
conflicting_bids_curr   =   cell(NUM_OF_PARTICLES,1);
full_allocations_explored_no_th = num_of_bids;
for j = 2 : NUM_OF_PARTICLES
    full_allocations_explored_no_th = full_allocations_explored_no_th * (num_of_bids - PARTICLE_WIDTH_x * PARTICLE_WIDTH_y); 
end
full_allocations_explored_no_th = full_allocations_explored_no_th/factorial(NUM_OF_PARTICLES);
full_allocations_explored_with_th = 0;
num_of_pruning_conditions_calc = 0;
first_allocation = 1;
for i = 1 : num_of_bids
    if NUM_OF_PARTICLES == 1
        opt_revenue = max(sorted_bids_array(:,2));
        opt_comb_allocation = sorted_bids_array(sorted_bids_array(:,2) == opt_revenue,1);
        naive_allocation = opt_comb_allocation;
        naive_revenue = opt_revenue;
        full_allocations_explored_with_th = 0;
        num_of_pruning_conditions_calc = 0;
        break;
    end
    num_of_bids_to_explore  =   ones(NUM_OF_PARTICLES-1,1);
    count                   =   zeros(NUM_OF_PARTICLES-1,1);
    has_visited_level       =   zeros(NUM_OF_PARTICLES,1);
    curr_bid_num = i;
    level = 1;
    curr_alloc = sorted_bids_array(curr_bid_num,:);
    
    conflicting_bids_curr{level} = conflicting_bids_with_curr_alloc_modified(curr_alloc,conflicting_bids_array,PARTICLE_WIDTH_x,PARTICLE_WIDTH_y,[]);
    if length(conflicting_bids_curr{level}) >= length(sorted_bids_array(:,2))
        continue;
    end
    bids_array_no_i{level} = remove_conflicting_bids(sorted_bids_array,conflicting_bids_curr{level});
    bids_array_no_i{level} = bids_array_no_i{level}(~isnan(bids_array_no_i{level}(:,1)),:);
    bids_left_num{level} = bids_array_no_i{level}(:,1);
    if ~has_visited_level(level) 
        price_curr_alloc = sum(curr_alloc(:,2));
        num_of_bids_to_explore(level) = length(bids_left_num{level});
        has_visited_level(level) = 1;
        h_pai = sum(bids_array_no_i{level}(1:NUM_OF_PARTICLES-1,2));
        num_of_pruning_conditions_calc = num_of_pruning_conditions_calc + 1;
        if price_curr_alloc + h_pai < opt_revenue
            continue;
        end
    end
    %%
    one_level_up_flag = 0;
    while num_of_bids_to_explore(1)
        if one_level_up_flag
            if isempty(bids_left_num{1})
                break;
            end
            curr_alloc = curr_alloc(1:end-1,:);
            bids_left_num{level} = bids_left_num{level}(2:end);
            count(level) = count(level) + 1;
            one_level_up_flag = 0;
            if isempty(bids_left_num{length(curr_alloc(:,1))})
                curr_alloc = curr_alloc(1:end-1,:);
                if isempty(curr_alloc)
                    break;
                end
                bids_left_num{level - 1} = bids_left_num{level - 1}(2:end);
                count(level:end) = 0;
                has_visited_level(level:end) = 0;
                level = level - 1;
                count(level) = count(level) + 1;
            end
        end        
        x = ~(count - num_of_bids_to_explore);
        if sum(x) %one of them is zero
            %We've explored all bids, need to backtrack one level up
            fprintf("hit")
            one_level_up_flag = 1;
            has_visited_level(level) = 0;
            count(level) = 0;
            if level == 1
                break
            end
        end
        if (isempty(bids_left_num{level}))||(one_level_up_flag && isempty(bids_left_num{level-1})) %%WORKED so far without the or condition 
            curr_alloc = curr_alloc(1:end-1,:);
            if isempty(curr_alloc)
                break;
            end
            bids_left_num{level - 1} = bids_left_num{level - 1}(2:end);
            count(level:end) = 0;
            has_visited_level(level:end) = 0;
            level = level - 1;
            count(level) = count(level) + 1;
        end

        valid_bid_idx = bids_left_num{level}(1);

        while length(curr_alloc(:,1)) ~= NUM_OF_PARTICLES
            mask = bids_array_no_i{level}(bids_array_no_i{level}(:,1) == valid_bid_idx,:);
            curr_alloc = [curr_alloc ; mask];
            if length(curr_alloc(:,1)) == NUM_OF_PARTICLES
                break;
            end
            level = length(curr_alloc(:,1));
            conflicting_bids_curr{level} = conflicting_bids_with_curr_alloc_modified(curr_alloc,conflicting_bids_array,PARTICLE_WIDTH_x,PARTICLE_WIDTH_y,conflicting_bids_curr{level-1});
            bids_array_no_i{level} = remove_conflicting_bids(sorted_bids_array,conflicting_bids_curr{level});
            bids_array_no_i{level} = bids_array_no_i{level}(~isnan(bids_array_no_i{level}(:,1)),:);
            bids_left_num{level} = bids_array_no_i{level}(:,1);

            if ~has_visited_level(level)
                price_curr_alloc = sum(curr_alloc(:,2));
                num_of_bids_to_explore(level) = length(bids_left_num{level});
                has_visited_level(level) = 1;
                h_pai = sum(bids_array_no_i{level}(1:min(NUM_OF_PARTICLES-level,num_of_bids_to_explore(level)),2));
                num_of_pruning_conditions_calc = num_of_pruning_conditions_calc + 1;
                if price_curr_alloc + h_pai <= opt_revenue
                    bids_left_num{level} = [];
                    one_level_up_flag = 1;
                    has_visited_level(level) = 0;
                    count(level) = 0;
                    level = level - 1;
                end
            end
            if ~isempty(bids_left_num{length(curr_alloc(:,1))})
                valid_bid_idx = bids_left_num{length(curr_alloc(:,1))}(1);
            else
                break;
            end
        end
        if length(curr_alloc(:,1)) == NUM_OF_PARTICLES
            % need to explore all options for the last bid, since we didn't
            % prune before.
            % valid_bid_idx holds all bids that need to be explored
            full_allocations_explored_with_th = full_allocations_explored_with_th + 1;
            count(level) = count(level) + 1;
            iter_revenue = sum(curr_alloc(:,2));
            if iter_revenue > opt_revenue
                opt_revenue = iter_revenue;
                opt_comb_allocation = curr_alloc(:,1);
            end
            if first_allocation
                naive_allocation = curr_alloc(:,1);
                naive_revenue = iter_revenue;
                first_allocation = 0;
            end
            curr_alloc = curr_alloc(1:end-1,:); %remove last bid
            bids_left_num{level} = bids_left_num{level}(2:end);
            h_pai = sum(bids_array_no_i{level}(1:NUM_OF_PARTICLES-level,2));
            num_of_pruning_conditions_calc = num_of_pruning_conditions_calc + 1;
                if price_curr_alloc + h_pai <= opt_revenue
                    bids_left_num{level} = [];
                    one_level_up_flag = 1;
                    has_visited_level(level) = 0;
                    count(level) = 0;
                    level = level - 1;
                end
        end
    end


end


run_time = toc;
% fprintf("With row sorting the bids by maximal price: \n")
% fprintf("    num of allocations that should be explored without th: %d \n",full_allocations_explored_no_th);
% fprintf("    num of allocations explored with th: %d \n",full_allocations_explored_with_th);
% fprintf("    magnitude of improvement: %g \n",full_allocations_explored_no_th/full_allocations_explored_with_th);
% fprintf("    run time is : %g seconds \n",run_time);
% fprintf("    revenue is : %g \n",opt_revenue);

end