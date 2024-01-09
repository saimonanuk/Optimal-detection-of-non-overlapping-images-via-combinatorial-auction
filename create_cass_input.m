function [auction_sum,new_xcorr_pic_mat,bids_array] = create_cass_input(Particle_upper_left_corner_first_pic,xcorr_pic,PARTICLE_WIDTH_x,PARTICLE_WIDTH_y,orig_pic_size_rows,orig_pic_size_cols,cryo_test)
 if cryo_test
    new_xcorr_pic_mat = xcorr_pic;
 else
     new_xcorr_pic_mat = xcorr_pic(PARTICLE_WIDTH_y:end, PARTICLE_WIDTH_x:end); %The first particle we would possibly allocate is at [particle_width-1,particle_width-1]
 end

goods_array = zeros(orig_pic_size_rows * orig_pic_size_cols, PARTICLE_WIDTH_x * PARTICLE_WIDTH_y); %square particle, each pixel can be a particle left corner

p_ij_vec = reshape(new_xcorr_pic_mat,[],1).'; %(i,j) element in xcorr_pic -> (orig_pic_size_rows)*(j-1) + i

for p = 1 : length(p_ij_vec)
    orig_i = mod(p,orig_pic_size_rows);
    orig_j = floor(p/orig_pic_size_rows) + 1;
    k = 1;
    for q = 1 : PARTICLE_WIDTH_x
        goods_array(p, k:k + PARTICLE_WIDTH_y - 1) = orig_pic_size_rows * (q-1) + (orig_j-1)*orig_pic_size_rows + orig_i : orig_pic_size_rows * (q-1) + (orig_j-1)*orig_pic_size_rows + orig_i + PARTICLE_WIDTH_y - 1;
        k = k + PARTICLE_WIDTH_y;
    end
end

goods_array(goods_array > orig_pic_size_rows * orig_pic_size_cols) = NaN;
auction_sum = 0;

if ~cryo_test
    p_ij_vec = floor(10^4 * p_ij_vec); %%prices are integeres
end
num_of_goods = orig_pic_size_rows * orig_pic_size_cols; %num of pixels in the picture
num_of_bids = orig_pic_size_rows * orig_pic_size_cols; %every pixel can be estimated as a particle upper left corner


bids_vec = (1:num_of_bids).';
bids_array = [bids_vec,p_ij_vec.',goods_array];

