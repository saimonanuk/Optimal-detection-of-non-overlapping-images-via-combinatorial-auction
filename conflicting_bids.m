function [conflicting_bids_array] = conflicting_bids(bids_array,PARTICLE_WIDTH_x,PARTICLE_WIDTH_y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% - INPUT
%%%%%%          bids_array - an array of all bids and the goods requested
%%%%%%          by those bids
%%%%%%          Particle_width - the width of the particle in the picture -
%%%%%%          determines maximal goods number requested by each bid
%%%%%% - OUTPUT
%%%%%%          conflicting_bids_array - an array that is built such that:
%%%%%%          conflicting_bids_array(:,1) - good idx's in picture
%%%%%%          conflicting_bids_array(:,2:end) - the bid num that this
%%%%%%          good is requested by.
%%%%%%
%%%%%%  for example = all the bids wanting good number 7, are located in
%%%%%%  conflicting_bids_array(7,2:end).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_of_bids = length(bids_array(:,1)); %num_of_goods = num_of_bids
num_of_goods = num_of_bids;
%  This function calculates all conflicted bids
bids_array = [bids_array(:,1),bids_array(:,3:end)]; %The revenue of each bid is not interesting in the calculation of conflicting bids
bids_array_without_nan = fillmissing(bids_array,"constant",-1);
conflicting_bids_array = NaN(num_of_goods, 1 + PARTICLE_WIDTH_x * PARTICLE_WIDTH_y);
conflicting_bids_array(:,1) = (1:num_of_goods).'; %good_num
for i = 1 : num_of_goods
    a = find(~(bids_array_without_nan-i));
    a = a(2:end);
    a = mod(a,num_of_goods);
    conflicting_bids_array(i,2: 1+length(a)) = a;
end
conflicting_bids_array(num_of_goods,2) = num_of_goods;
%We make use of the fact we have num_of_goods = num_of_bids
end