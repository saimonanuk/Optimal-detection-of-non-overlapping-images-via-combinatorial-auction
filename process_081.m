clear; close all; clc;
map = load("empiar_10081_pic.mat");
map = map.ctf_corrected_micrograph;
avg_map = mean(map,'all');
new_map = map - avg_map;
figure;

noise_patch_1 = new_map(2250:2500-1,3300:3300+265-1); %n = 750 samples, p = 256
noise_patch_2 = new_map(1650:1650-1,3300:3300+265-1);
noise_patch_3 = new_map(2000:2400,1150:1150+265-1);
noise_patch_4 = new_map(1950:2200,444:444+265-1);
noise_patch_5 = new_map(1950:2200,720:720+265-1);
total_noise_patch = [noise_patch_1;noise_patch_2;noise_patch_3;noise_patch_4;noise_patch_5];
num_noise_samples = length(total_noise_patch(:,1));
noise_vec_size = length(noise_patch_1(1,:));
noise_mean_vec = mean(total_noise_patch,1).';
cov_mat = zeros(noise_vec_size,noise_vec_size);
for i = 1 : num_noise_samples
    cov_mat = cov_mat + (total_noise_patch(i,:).' - noise_mean_vec) * (total_noise_patch(i,:).' - noise_mean_vec).';
end
cov_mat = cov_mat / (num_noise_samples-1);
new_map_whitened = zeros(3710,3710);
cov_mat_inv = sqrtm(inv(cov_mat));
for i = 1 :length(new_map(:,1))
    for j = 1 : length(new_map(1,:))/265
        a = cov_mat_inv * new_map(i,(j-1)*265 + 1 : j*265).';
        new_map_whitened(i,(j-1)*265 + 1 : j*265) = a.';
    end
end
new_map_whitened = 10*new_map_whitened;
load dec_filter_10081_2d.mat
whitened_image_diff_way = imfilter(new_map_whitened,dec_filter_10081_2d,'same');
map_before_decimation_whitened = whitened_image_diff_way;
imagesc(map_before_decimation_whitened)
colormap("gray")
imcontrast()
decimated_image_with_whitening = map_before_decimation_whitened(1:12:end,1:12:end);
imagesc(decimated_image_with_whitening)
colormap("gray")
imcontrast
radius = 3;

se = strel('disk', radius);
binaryImage = zeros(2*radius+1, 2*radius+1);
[x, y] = meshgrid(1:size(binaryImage, 2), 1:size(binaryImage, 1));
binaryImage(sqrt((x-radius-1).^2 + (y-radius-1).^2) <= radius) = 1;
Particle = binaryImage;
total_pic_corr = (xcorr2(decimated_image_with_whitening,Particle));
figure;
imagesc(total_pic_corr)
colormap('gray')
imcontrast

image = decimated_image_with_whitening;
image(187:198,12:24) = 0;
image_corr = total_pic_corr;
image_corr(186:202,13:27) = 0;
figure;
imagesc(image_corr)
colormap('gray')
imcontrast
edge_map_corr = edge(image_corr,'prewitt');
figure;
imagesc(edge_map_corr)
colormap('gray')
new_edge_map = edge(image,'prewitt');
figure;
imagesc(new_edge_map)
colormap('gray')
edgeMap = edge_map_corr;

se = strel('disk', 3); 
closeEdgesMap = imclose(edgeMap, se);

labeledParticles = bwlabel(closeEdgesMap);
 
minParticleSize = 50; 
filteredParticles = labeledParticles;

stats = regionprops(filteredParticles, 'Area');
for i = 1:numel(stats)
    if stats(i).Area < minParticleSize
        filteredParticles(filteredParticles == i) = 0;
    end
end
% Display
figure;
subplot(1, 2, 1);
imagesc(new_edge_map);
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
crossCorrelation = total_pic_corr;
threshold = 48;
binaryMask = crossCorrelation > threshold; 
specificBox = [257 - 1, 110 - 1, 9 + 2, 10 + 2]; 

x = round(specificBox(1));
y = round(specificBox(2));
width = round(specificBox(3));
height = round(specificBox(4));
roiBinaryMask = binaryMask(y:y+height, x:x+width);
particleCount = max(max(bwlabel(roiBinaryMask)));
labeledParticles = bwlabel(roiBinaryMask);
stats = regionprops(labeledParticles, 'BoundingBox');
particleCoordinates = cat(1, stats.BoundingBox);
newImage = roiBinaryMask; 
particleRadius = 3;
for i = 1:size(particleCoordinates, 1)
    xLeftUpper = round(particleCoordinates(i, 1));
    yLeftUpper = round(particleCoordinates(i, 2));
    if xLeftUpper >= particleRadius && yLeftUpper >= particleRadius && ...
       (xLeftUpper + 2*particleRadius) <= size(newImage, 2) && ...
       (yLeftUpper + 2*particleRadius) <= size(newImage, 1)
        
        newImage = insertShape(double(newImage), 'FilledCircle', ...
            [xLeftUpper + particleRadius, yLeftUpper + particleRadius, particleRadius], 'Color', 'red');
    end
end
test_image = image(y-3:y+height-3, x-3:x+width-3);
pic_corr = total_pic_corr(y:y+height, x:x+width);
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
%%
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

estimated_image_naive = place_particles_left_corner(particle_estimation_naive - [3 3], PICTURE_SIZE_ROWS, PICTURE_SIZE_COLS, PARTICLE_WIDTH_x, PARTICLE_WIDTH_y, Particle);
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
estimated_image_comb = place_particles_left_corner(particle_estimation_comb - [3 3], PICTURE_SIZE_ROWS, PICTURE_SIZE_COLS, PARTICLE_WIDTH_x, PARTICLE_WIDTH_y, Particle);
estimated_image_comb_figure = estimated_image_comb;
estimated_image_comb_figure(estimated_image_comb_figure > 0) = 1;
figure;
imagesc(imfilter(estimated_image_comb_figure,ones(2,2),'same'))
colormap('gray')
imcontrast
figure;
imagesc(imfilter(test_image,ones(3,3),'same'))
colormap('gray')
imcontrast