function [Particle_centers] = generate_particles_upper_left_corner_far(Num_of_particles,PARTICLE_WIDTH_x,PARTICLE_WIDTH_y,orig_pic_size_rows,orig_pic_size_cols, DISTANCE)

%%%% The particle centers generated in this function will be DISTANCE apart
% from each other. while Distance >= 2*(Particle_width) - 1 We can
% distinguish between the particles perfectly. as Distance approches
% Distance = Particle_width we will experience worse performance.
rng('shuffle')

Particle_centers = [];
size_particle_centers = 0;
max_particle_width = max(PARTICLE_WIDTH_x,PARTICLE_WIDTH_y);

if DISTANCE >= (2*(max_particle_width) - 1) %% Generate particles, we can distinguish between them
    while size_particle_centers < Num_of_particles
        new_center = [randi(orig_pic_size_rows - PARTICLE_WIDTH_y + 1, 1 , 1) randi(orig_pic_size_cols - PARTICLE_WIDTH_x + 1, 1 , 1)]; %Whole particle has to be in the picture
        if isempty(Particle_centers)
            Particle_centers = new_center;
        else
            norm_closest_center = vecnorm(new_center - Particle_centers,2,2);
            if ~sum(~(norm_closest_center >= DISTANCE))
                Particle_centers = [Particle_centers ; new_center];
            end
        end
        size_particle_centers = size(Particle_centers,1);
    end
elseif  (DISTANCE >= max_particle_width + 1) && (DISTANCE < ((2*(max_particle_width) - 1)))
    while size_particle_centers(1) < Num_of_particles
        new_center = [randi(orig_pic_size_rows - PARTICLE_WIDTH_y + 1, 1 , 1) randi(orig_pic_size_cols - PARTICLE_WIDTH_x + 1, 1 , 1)];
        if isempty(Particle_centers)
            new_center = [randi(20 - PARTICLE_WIDTH_y + 1, 1 , 1) randi(20 - PARTICLE_WIDTH_x + 1, 1 , 1)];
            Particle_centers = new_center;
        else
            norm_closest_center = vecnorm(new_center - Particle_centers,2,2);
            if ~sum(~((norm_closest_center >= DISTANCE) & (norm_closest_center <= DISTANCE + 1)))
                Particle_centers = [Particle_centers ; new_center];
            end
        end
        size_particle_centers = size(Particle_centers,1);
    end
else
    %%%%%%%%% Distance is exactly particle width, particles aren't
    %%%%%%%%% seperated at all (but aren't overlapping)
    Particle_centers = [randi(orig_pic_size_rows- PARTICLE_WIDTH_y + 1, 1 , 1) randi(orig_pic_size_cols- PARTICLE_WIDTH_x + 1, 1 , 1)];

    size_particle_centers = size(Particle_centers,1);
    right_left_up_down_array = [0,0]; 
 
    if Particle_centers(end,2) < orig_pic_size_cols - 2*PARTICLE_WIDTH_x
        right_left_up_down_array(1) = 1;
    else
        right_left_up_down_array(1) = 0;
    end
    if Particle_centers(end,1) >  PARTICLE_WIDTH_y 
        right_left_up_down_array(2) = 1; % go up
    else
        right_left_up_down_array(2) = 0; % go down
    end

    while size_particle_centers < Num_of_particles
        if right_left_up_down_array(1)
           new_center = [Particle_centers(end,1), Particle_centers(end,2) +  PARTICLE_WIDTH_x];
        else
           new_center = [Particle_centers(end,1), Particle_centers(end,2) -  PARTICLE_WIDTH_x];
        end
        Particle_centers = [Particle_centers ; new_center];
        size_particle_centers = size(Particle_centers,1);
        if size_particle_centers == Num_of_particles
            break;
        end
        if right_left_up_down_array(1) && Particle_centers(end,2) >= orig_pic_size_cols - 2*PARTICLE_WIDTH_x
            %% need to move up or down based on right_left_up_down_array(2), and then start to move left
            new_center_row_shift = sign(2*right_left_up_down_array(2) - 1) * PARTICLE_WIDTH_y;
            new_center = [Particle_centers(end,1) - new_center_row_shift, Particle_centers(end,2)];
            Particle_centers = [Particle_centers ; new_center];
            size_particle_centers = size(Particle_centers,1);
            right_left_up_down_array(1) = 0;
        elseif  ~right_left_up_down_array(1) && Particle_centers(end,2) <= PARTICLE_WIDTH_x
            new_center_row_shift = sign(2*right_left_up_down_array(2) - 1) * PARTICLE_WIDTH_y;
            new_center = [Particle_centers(end,1) - new_center_row_shift, Particle_centers(end,2)];
            Particle_centers = [Particle_centers ; new_center];
            size_particle_centers = size(Particle_centers,1);
            right_left_up_down_array(1) = 1;
        end
    end
end
Particle_centers = sortrows(Particle_centers); %In large particle widths, not always NUM_OF_PARTICLES particles are generated. problem?
end