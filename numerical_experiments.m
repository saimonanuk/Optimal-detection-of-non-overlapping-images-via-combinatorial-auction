clear; clc; close all;
NUM_OF_RUNS             = 1000;
PICTURE_SIZE_ROWS       = 40;
PICTURE_SIZE_COLS       = 40;
Num_of_particles_in_pic = (4).';
PARTICLE_WIDTH              = 3;
PARTICLE_WIDTH_x            = 3;
PARTICLE_WIDTH_y            = 3;
Particle_distance       = PARTICLE_WIDTH;
SNR                     =   (8:-2:0).';
N_0                     = 10.^(-SNR/10);
Particle                = ones(PARTICLE_WIDTH_y,PARTICLE_WIDTH_x);
Precision_comb = 0;
Recall_comb = 0;
F_1_comb = 0;
Precision_naive = 0;
Recall_naive = 0;
F_1_naive = 0;
full_allocations_explored = zeros(length(SNR),1);
num_of_prun_cond_calc = zeros(length(SNR),1);
run_time = zeros(length(SNR),1);
TP_comb = zeros(length(SNR),1);
TP_naive = zeros(length(SNR),1);
for p = 1 : length(SNR)
    for k = 1 : NUM_OF_RUNS
        gaussian_noise = sqrt(N_0(p))*(randn(PICTURE_SIZE_ROWS,PICTURE_SIZE_COLS));
        picture = gaussian_noise;
        Particle_upper_left_corner_first_pic  = generate_particles_upper_left_corner_far(Num_of_particles_in_pic,PARTICLE_WIDTH_x,PARTICLE_WIDTH_y,PICTURE_SIZE_ROWS,PICTURE_SIZE_COLS, Particle_distance);
        first_pic_with_particles    = place_particles_left_corner(Particle_upper_left_corner_first_pic, PICTURE_SIZE_ROWS, PICTURE_SIZE_COLS, PARTICLE_WIDTH_x, PARTICLE_WIDTH_y, Particle);
        picture = picture + first_pic_with_particles;
        xcorr_pic_Particle = xcorr2(picture, Particle);
        [auction_sum, new_xcorr_mat,bids_array] = create_cass_input(Particle_upper_left_corner_first_pic,xcorr_pic_Particle,PARTICLE_WIDTH_x,PARTICLE_WIDTH_y,PICTURE_SIZE_ROWS,PICTURE_SIZE_COLS,0);
        [opt_alloc_with_sort, opt_rev_with_sort,run_time_with_sort,naive_allocation,naive_revenue,full_allocations_explored_with_th,num_of_prun_cond_calc_iter] = find_opt_allocation_sorted_bids_modified(bids_array,Num_of_particles_in_pic,PARTICLE_WIDTH_x,PARTICLE_WIDTH_y);
        winning_bids_vec_naive = naive_allocation;
        particle_estimation_naive = zeros(length(winning_bids_vec_naive),2);
        for q = 1 : length(winning_bids_vec_naive)
            particle_estimation_naive(q,:) = [mod(winning_bids_vec_naive(q),PICTURE_SIZE_ROWS) floor(winning_bids_vec_naive(q)/PICTURE_SIZE_ROWS) + 1] - [7 7];
        end
        particle_estimation_naive = sortrows(particle_estimation_naive);
        particle_estimation_naive = particle_estimation_naive + [2*PARTICLE_WIDTH_y + 1, 2*PARTICLE_WIDTH_x + 1];
        winning_bids_vec_comb = opt_alloc_with_sort;
        particle_estimation_comb = zeros(length(winning_bids_vec_comb),2);
        for q = 1 : length(winning_bids_vec_comb)
            particle_estimation_comb(q,:) = [mod(winning_bids_vec_comb(q),PICTURE_SIZE_ROWS) floor(winning_bids_vec_comb(q)/PICTURE_SIZE_ROWS) + 1] -[7 7];
        end
        particle_estimation_comb = sortrows(particle_estimation_comb);
        particle_estimation_comb = particle_estimation_comb + [2*PARTICLE_WIDTH_y + 1, 2*PARTICLE_WIDTH_x + 1];
        TP_comb_iter = TP_calc(Particle_upper_left_corner_first_pic,particle_estimation_comb,Num_of_particles_in_pic,PARTICLE_WIDTH_y,PARTICLE_WIDTH_x);
        TP_naive_iter = TP_calc(Particle_upper_left_corner_first_pic,particle_estimation_naive,Num_of_particles_in_pic,PARTICLE_WIDTH_y,PARTICLE_WIDTH_x);
        TP_comb(p) = TP_comb(p) + TP_comb_iter;
        TP_naive(p) = TP_naive(p) + TP_naive_iter;
        % Print the current status
        fprintf('Running for p = %d, k = %d \n', p, k);
    end
end
    TP_comb = TP_comb/NUM_OF_RUNS;
    TP_naive = TP_naive/NUM_OF_RUNS;





