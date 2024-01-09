clear; clc; close all;
map = ReadMRC("18apr21a_2Bsr_00005gr_00005sq_v02_00002hln_00002esn-a-DW_ctf_pf_10217.mrcs");
avg_map = mean(map,'all');
new_map = map - avg_map;
load dec_filter_10081_2d.mat
image_after_filt = imfilter(100*new_map,dec_filter_10081_2d,'same');
imagesc(image_after_filt)
colormap('gray')
imcontrast
figure;
decimated_image = image_after_filt(1:12:end,1:12:end)/10;
imagesc(decimated_image)
colormap('gray')
imcontrast
figure;
avg_filt_image = imfilter(decimated_image,ones(2,2),'same');
imagesc(avg_filt_image)
colormap('gray')
imcontrast


radius = 6;
se = strel('disk', radius);
binaryImage = zeros(2*radius+1, 2*radius+1);
[x, y] = meshgrid(1:size(binaryImage, 2), 1:size(binaryImage, 1));
binaryImage(sqrt((x-radius-1).^2 + (y-radius-1).^2) <= radius) = 1;
Particle_6 = binaryImage;
imcontrast
estimated_particle1 = decimated_image(35:46,162:172);
estimated_particle2 = decimated_image(111:122,228:238);
estimated_particle3 = decimated_image(41:52,210:220);
estimated_particle4 = decimated_image(23:34,170:180);

estimated_particle = estimated_particle1 + estimated_particle2 + estimated_particle3 + estimated_particle4;
total_pic_corr_no_abs = (xcorr2(decimated_image,Particle_6));
figure;
imagesc(total_pic_corr_no_abs)
colormap('gray')
imcontrast
image = decimated_image;
image_corr = total_pic_corr_no_abs;
edge_map_corr = edge(image_corr,'prewitt');
figure;
imagesc(edge_map_corr)
colormap('gray')
closeEdgesMap = imclose(edge_map_corr, se);
labeledParticles = bwlabel(closeEdgesMap);
minParticleSize = 120; 
filteredParticles = labeledParticles;
stats = regionprops(filteredParticles, 'Area');
for i = 1:numel(stats)
    if stats(i).Area < minParticleSize
        filteredParticles(filteredParticles == i) = 0;
    end
end
figure;
subplot(1, 2, 1);
imagesc(edge_map_corr);
title('Original Image');
subplot(1, 2, 2);
imagesc(filteredParticles);
title('Filtered Image');
colormap('gray')
stats = regionprops(filteredParticles, 'BoundingBox');

closeParticlesRegionImage = zeros(size(filteredParticles));
for i = 1:numel(stats)
    bbox = round(stats(i).BoundingBox);
    xmin = max(1, bbox(1)); 
    ymin = max(1, bbox(2));
    xmax = min(size(filteredParticles, 2), xmin + bbox(3));
    ymax = min(size(filteredParticles, 1), ymin + bbox(4));

    closeParticlesRegionImage(ymin:ymax, xmin:xmax) = 1;
end
figure;
imagesc(closeParticlesRegionImage);
colormap('gray')

crossCorrelation = image_corr;
threshold = 18; 
binaryMask = crossCorrelation > threshold;
specificBox = [50 - 3, 215 - 3, 12 + 6, 19 + 6];
x = round(specificBox(1));
y = round(specificBox(2));
width = round(specificBox(3));
height = round(specificBox(4));

roiBinaryMask = binaryMask(y:y+height, x:x+width);
particleCount = max(max(bwlabel(roiBinaryMask)));
labeledParticles = bwlabel(roiBinaryMask);


test_image = image(y-6:y+height-6, x-6:x+width-6);
pic_corr = image_corr(y:y+height, x:x+width);
Particle = Particle_6;

xcorr_pic_Particle = (pic_corr);
PARTICLE_WIDTH_x = length(Particle(1,:));
PARTICLE_WIDTH_y = length(Particle(:,1));
PICTURE_SIZE_COLS = length(test_image(1,:));
PICTURE_SIZE_ROWS = length(test_image(:,1));
image_area = PICTURE_SIZE_COLS * PICTURE_SIZE_ROWS;
particle_area = PARTICLE_WIDTH_x * PARTICLE_WIDTH_y;
max_num_of_particles_in_image = floor(image_area/particle_area - 1/2); 
k_values = 1:max_num_of_particles_in_image;
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
a = 0;
for k = 1 : length(k_values)
    [new_xcorr_mat,bids_array] = create_cass_input_cryo(xcorr_pic_Particle,PARTICLE_WIDTH_x,PARTICLE_WIDTH_y,PICTURE_SIZE_ROWS,PICTURE_SIZE_COLS);
    [opt_comb_allocation{k}, opt_revenue(k),run_time(k),naive_allocation{k}, naive_revenue(k), full_allocations_explored_with_th, num_of_pruning_conditions_calc] = find_opt_allocation_sorted_bids_modified(bids_array,k_values(k),PARTICLE_WIDTH_x,PARTICLE_WIDTH_y);
    iterations_made = P;
    for i = 1 : P
        linear_test_image = test_image(:);
        curr_permute = randperm(length(linear_test_image)).';
        linear_test_image_P = linear_test_image(curr_permute);
        curr_test_image_mat = reshape(linear_test_image_P,[PICTURE_SIZE_ROWS,PICTURE_SIZE_COLS]);
        curr_xcorr_mat = (xcorr2(curr_test_image_mat,Particle));
        [~,bids_array_p] = create_cass_input_cryo(curr_xcorr_mat(PARTICLE_WIDTH_y:end,PARTICLE_WIDTH_x:end),PARTICLE_WIDTH_x,PARTICLE_WIDTH_y,PICTURE_SIZE_ROWS,PICTURE_SIZE_COLS);
        [~, opt_revenue_p, run_time_p, ~ , naive_revenue_p , ~, ~] = find_opt_allocation_sorted_bids_modified(bids_array_p,k_values(k),PARTICLE_WIDTH_x,PARTICLE_WIDTH_y);
        opt_rev_noise(k) = opt_rev_noise(k) + opt_revenue_p;
        naive_rev_noise(k) = naive_rev_noise(k) + naive_revenue_p;
        run_time_noise(k) = run_time_noise(k) + run_time_p;
    end
    avg_run_time(k) = run_time_noise(k)/iterations_made;
    opt_gap_statistics(k) = -opt_rev_noise(k)/iterations_made + opt_revenue(k);
    naive_gap_statistics(k) = -naive_rev_noise(k)/iterations_made + naive_revenue(k);
end

[max_gap,k_opt] = max(abs(opt_gap_statistics));
[max_gap_naive,k_opt_naive] = max(abs(naive_gap_statistics));
opt_comb_allocation_total = opt_comb_allocation{k_opt};
opt_revenue_total = opt_revenue(k_opt);
naive_allocation_total = naive_allocation{k_opt_naive};
naive_revenue_total = naive_revenue(k_opt_naive);


winning_bids_vec_naive = naive_allocation_total;
particle_estimation_naive = zeros(length(winning_bids_vec_naive),2);
for q = 1 : length(winning_bids_vec_naive)
    particle_estimation_naive(q,:) = [mod(winning_bids_vec_naive(q),PICTURE_SIZE_ROWS) floor(winning_bids_vec_naive(q)/PICTURE_SIZE_ROWS) + 1];
end

estimated_image_naive = place_particles_left_corner(particle_estimation_naive - [6  6], PICTURE_SIZE_ROWS, PICTURE_SIZE_COLS, PARTICLE_WIDTH_x, PARTICLE_WIDTH_y, Particle);
estimated_image_naive_figure = estimated_image_naive;
estimated_image_naive_figure(estimated_image_naive_figure > 0) = 1;
figure;
imagesc(imfilter(estimated_image_naive_figure,ones(2,2),'same'))
colormap('gray')
imcontrast

winning_bids_vec_comb = opt_comb_allocation_total;
particle_estimation_comb = zeros(length(winning_bids_vec_comb),2);
           
for q = 1 : length(winning_bids_vec_comb)
    particle_estimation_comb(q,:) = [mod(winning_bids_vec_comb(q),PICTURE_SIZE_ROWS) floor(winning_bids_vec_comb(q)/PICTURE_SIZE_ROWS) + 1];
end
particle_estimation_comb = sortrows(particle_estimation_comb);
estimated_image_comb = place_particles_left_corner(particle_estimation_comb - [6 6], PICTURE_SIZE_ROWS, PICTURE_SIZE_COLS, PARTICLE_WIDTH_x, PARTICLE_WIDTH_y, Particle);
estimated_image_comb_figure = estimated_image_comb;
estimated_image_comb_figure(estimated_image_comb_figure > 0) = 1;
figure;
imagesc(imfilter(estimated_image_comb_figure,ones(2,2),'same'))
colormap('gray')
imcontrast
figure;
imagesc(imfilter(test_image,ones(5,5),'same'))
colormap('gray')
imcontrast
