function [picture] = place_particles_left_corner(Particle_centers,Picture_size_rows, Picture_size_cols,PARTICLE_WIDTH_x,PARTICLE_WIDTH_y, Particle)
            picture = zeros(Picture_size_rows,Picture_size_cols);
        for i = 1:length(Particle_centers(:,1)) % generate particle in picture for each center
        
            Particle_y_loc = (Particle_centers(i,1): Particle_centers(i,1) + PARTICLE_WIDTH_y - 1);
            Particle_x_loc = (Particle_centers(i,2): Particle_centers(i,2) + PARTICLE_WIDTH_x - 1);
            Particle_valid_cond_x = (Particle_y_loc > 0) & (Particle_y_loc <= Picture_size_rows);
            Particle_valid_cond_y = (Particle_x_loc > 0) & (Particle_x_loc <= Picture_size_cols);
            valid_Particle_y_loc = Particle_y_loc(Particle_valid_cond_x); % valid x idx's
            valid_Particle_x_loc = Particle_x_loc(Particle_valid_cond_y); % valid y idx's
            Particle_i = Particle(Particle_valid_cond_x,Particle_valid_cond_y); % add Partial particle to the picture
            picture(valid_Particle_y_loc,valid_Particle_x_loc) = picture(valid_Particle_y_loc,valid_Particle_x_loc)...
                                                                 + Particle_i;
        end

end