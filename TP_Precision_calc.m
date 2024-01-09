function [T_P,Precision] = TP_Precision_calc(real_loc_particles,est_loc_particles,NUM_OF_PARTICLES,Particle_width_y,Particle_width_x)
    T_P = 0;
    est_loc_particle_size = length(est_loc_particles(:,1));
    remaining_real_loc_particles = real_loc_particles;
    if est_loc_particle_size < length(real_loc_particles(:,1))
        bq = 1;
    end
    for q = 1 : NUM_OF_PARTICLES
        particle_real_q = real_loc_particles(q,:);
        diff_norm_vec = vecnorm(particle_real_q - est_loc_particles , 2 , 2);
        k = find(~diff_norm_vec);
        if k 
            new_diff_vec = vecnorm(particle_real_q - remaining_real_loc_particles , 2 , 2);
            j = find(~new_diff_vec);
            est_loc_particles(k,:) = [];
            % successfull estimation
            remaining_real_loc_particles(j,:) = [];
            T_P = T_P + 1;
%             est_loc_particles(q,:) = [];
        end
    end
    new_remaining_real_loc = remaining_real_loc_particles;
    for i = 1 :length(remaining_real_loc_particles(:,1))
        particle_real_i = remaining_real_loc_particles(i,:);
        diff_vec = particle_real_i - est_loc_particles;
        dist_metric = all(abs(diff_vec) <= max(Particle_width_y,Particle_width_x)/2, 2);
        exists = any(dist_metric);
        indices = find(dist_metric);
        if exists
            new_diff_vec = vecnorm(particle_real_i - new_remaining_real_loc , 2 , 2);
            j = find(~new_diff_vec);
            est_loc_particles(indices(1),:) = [];
            new_remaining_real_loc(j,:) = [];
            T_P = T_P + 1;
        end
    end
    Precision = T_P/est_loc_particle_size;
    T_P = T_P/length(real_loc_particles(:,1));
end