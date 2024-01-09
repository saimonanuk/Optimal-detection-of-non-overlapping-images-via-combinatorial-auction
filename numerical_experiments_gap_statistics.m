clear; clc; close all;
NUM_OF_RUNS             = 500;
PICTURE_SIZE_ROWS       = 40;
PICTURE_SIZE_COLS       = 40;
Num_of_particles_in_pic = 4;
PARTICLE_WIDTH              = 3;
PARTICLE_WIDTH_x            = 3;
PARTICLE_WIDTH_y            = 3;
Particle_distance       = PARTICLE_WIDTH;
SNR                     =   (8:-2:0).';
N_0                     = 10.^(-SNR/10);
Particle                = ones(PARTICLE_WIDTH_y,PARTICLE_WIDTH_x);
Recall_comb = 0;
F_1_comb = 0;
Recall_naive = 0;
F_1_naive = 0;
full_allocations_explored = zeros(length(SNR),1);
num_of_prun_cond_calc = zeros(length(SNR),1);
TP_comb = zeros(length(SNR),1);
TP_naive = zeros(length(SNR),1);
Precision_comb = zeros(length(SNR),1);
Precision_naive = zeros(length(SNR),1);

image_area = PICTURE_SIZE_COLS * PICTURE_SIZE_ROWS;
particle_area = PARTICLE_WIDTH_x * PARTICLE_WIDTH_y;
max_num_of_particles_in_image = floor(image_area/particle_area - 1/2); 
k_values = 1:(max_num_of_particles_in_image);

opt_comb_allocation = cell(length(k_values),1);
opt_revenue         = zeros(length(k_values),1);
run_time            = zeros(length(k_values),1);
naive_allocation    = cell(length(k_values),1);
naive_revenue       = zeros(length(k_values),1);

P = 50;

opt_rev_noise = zeros(length(k_values),1);
naive_rev_noise = zeros(length(k_values),1);
run_time_noise = zeros(length(k_values),1);
opt_gap_statistics = zeros(length(k_values),1);
naive_gap_statistics = zeros(length(k_values),1);
avg_run_time = zeros(length(k_values),1);
k_opt_comb = zeros(length(NUM_OF_RUNS),1);
k_opt_naive = zeros(length(NUM_OF_RUNS),1);



for l = 1 : length(SNR)
    for p = 1 : NUM_OF_RUNS
        a = 0;
        gaussian_noise = sqrt(N_0(l))*(randn(PICTURE_SIZE_ROWS,PICTURE_SIZE_COLS));
        picture = gaussian_noise;
        Particle_upper_left_corner_first_pic  = generate_particles_upper_left_corner_far(Num_of_particles_in_pic,PARTICLE_WIDTH_x,PARTICLE_WIDTH_y,PICTURE_SIZE_ROWS,PICTURE_SIZE_COLS, Particle_distance);
        first_pic_with_particles    = place_particles_left_corner(Particle_upper_left_corner_first_pic, PICTURE_SIZE_ROWS, PICTURE_SIZE_COLS, PARTICLE_WIDTH_x, PARTICLE_WIDTH_y, Particle);
        picture = picture + first_pic_with_particles;

        xcorr_pic_Particle = xcorr2(picture, Particle);
        opt_rev_noise = zeros(length(k_values),1);
        naive_rev_noise = zeros(length(k_values),1);
        run_time_noise = zeros(length(k_values),1);
        opt_gap_statistics = zeros(length(k_values),1);
        naive_gap_statistics = zeros(length(k_values),1);
        avg_run_time = zeros(length(k_values),1);
        for k = 1 : length(k_values)
            [auction_sum, new_xcorr_mat,bids_array] = create_cass_input(Particle_upper_left_corner_first_pic,xcorr_pic_Particle,PARTICLE_WIDTH_x,PARTICLE_WIDTH_y,PICTURE_SIZE_ROWS,PICTURE_SIZE_COLS,0);
            [opt_comb_allocation{k}, opt_revenue(k),run_time(k),naive_allocation{k}, naive_revenue(k), full_allocations_explored_with_th, num_of_pruning_conditions_calc] = find_opt_allocation_sorted_bids_modified(bids_array,k_values(k),PARTICLE_WIDTH_x,PARTICLE_WIDTH_y);
            iterations_made = P;
            for i = 1 : P
                linear_test_image = picture(:);
                curr_permute = randperm(length(linear_test_image)).';
                linear_test_image_P = linear_test_image(curr_permute);
                curr_test_image_mat = reshape(linear_test_image_P,[PICTURE_SIZE_ROWS,PICTURE_SIZE_COLS]);
                curr_xcorr_mat = (xcorr2(curr_test_image_mat,Particle));
                [~,~,bids_array_p] = create_cass_input(Particle_upper_left_corner_first_pic,curr_xcorr_mat,PARTICLE_WIDTH_x,PARTICLE_WIDTH_y,PICTURE_SIZE_ROWS,PICTURE_SIZE_COLS,0);
                [~, opt_revenue_p, run_time_p, ~ , naive_revenue_p , ~, ~] = find_opt_allocation_sorted_bids_modified(bids_array_p,k_values(k),PARTICLE_WIDTH_x,PARTICLE_WIDTH_y);
                opt_rev_noise(k) = opt_rev_noise(k) + opt_revenue_p;
                naive_rev_noise(k) = naive_rev_noise(k) + naive_revenue_p;
                run_time_noise(k) = run_time_noise(k) + run_time_p;
                if (mod(i,8) == 0 && k >= Num_of_particles_in_pic)
                    opt_gap_statistics(k) = -opt_rev_noise(k)/i + opt_revenue(k);
                    if opt_gap_statistics(k) <= opt_gap_statistics(k-1)
                        iterations_made = i;
                        a = 1;
                        break;
                    end
                end
            end
            avg_run_time(k) = run_time_noise(k)/iterations_made;
            opt_gap_statistics(k) = -opt_rev_noise(k)/iterations_made + opt_revenue(k);
            naive_gap_statistics(k) = -naive_rev_noise(k)/iterations_made + naive_revenue(k);
        end
        [max_gap,k_opt_comb(p)] = max(abs(opt_gap_statistics));
        [max_gap_naive,k_opt_naive(p)] = max(abs(naive_gap_statistics));
        opt_comb_allocation_total = opt_comb_allocation{k_opt_comb};
        opt_revenue_total = opt_revenue(k_opt_comb);
        naive_revenue_total = naive_revenue(k_opt_naive);
        winning_bids_vec_naive = naive_allocation{k_opt_naive};
        particle_estimation_naive = zeros(length(winning_bids_vec_naive),2);
        for q = 1 : length(winning_bids_vec_naive)
            particle_estimation_naive(q,:) = [mod(winning_bids_vec_naive(q),PICTURE_SIZE_ROWS) floor(winning_bids_vec_naive(q)/PICTURE_SIZE_ROWS) + 1] - [7 7];
        end
        particle_estimation_naive = sortrows(particle_estimation_naive);
        particle_estimation_naive = particle_estimation_naive + [2*PARTICLE_WIDTH_y + 1, 2*PARTICLE_WIDTH_x + 1];
        winning_bids_vec_comb = opt_comb_allocation{k_opt_comb};
        particle_estimation_comb = zeros(length(winning_bids_vec_comb),2);
        for q = 1 : length(winning_bids_vec_comb)
            particle_estimation_comb(q,:) = [mod(winning_bids_vec_comb(q),PICTURE_SIZE_ROWS) floor(winning_bids_vec_comb(q)/PICTURE_SIZE_ROWS) + 1] -[7 7];
        end
        particle_estimation_comb = sortrows(particle_estimation_comb);
        particle_estimation_comb = particle_estimation_comb + [2*PARTICLE_WIDTH_y + 1, 2*PARTICLE_WIDTH_x + 1];
        [TP_comb_iter,Precision_comb_iter] = TP_Precision_calc(Particle_upper_left_corner_first_pic,particle_estimation_comb,k_opt_comb,PARTICLE_WIDTH_y,PARTICLE_WIDTH_x);
        [TP_naive_iter,Precision_naive_iter] = TP_Precision_calc(Particle_upper_left_corner_first_pic,particle_estimation_naive,k_opt_naive,PARTICLE_WIDTH_y,PARTICLE_WIDTH_x);
        TP_comb(l) = TP_comb(l) + TP_comb_iter;
        TP_naive(l) = TP_naive(l) + TP_naive_iter;
        Precision_comb(l) = Precision_comb(l) + Precision_comb_iter;
        Precision_naive(l) = Precision_naive(l) + Precision_naive_iter;
        fprintf('Running for l = %d, p = %d \n', l, p);

    end
end
    TP_comb = TP_comb/NUM_OF_RUNS;
    TP_naive = TP_naive/NUM_OF_RUNS;
    Precision_comb = Precision_comb/NUM_OF_RUNS;
    Precision_naive = Precision_naive/NUM_OF_RUNS;
    save('TP_comb_gap','TP_comb');
    save('TP_naive_gap','TP_naive');
    save('Precision_comb_gap','Precision_comb');
    save('Precision_naive_gap','Precision_naive');

