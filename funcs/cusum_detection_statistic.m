%function [detect_stat_seq] = kernel_detection_statistic(pre_change_sample,post_change_sample,omega_B,Num_blk,varargin)
function [detect_stat_seq] = cusum_detection_statistic(pre_change_sample,post_change_sample,sgma, delta)



    

    [sample_size_ref,~] = size(pre_change_sample);
    [sample_size,~] = size(post_change_sample);
    
    
    
    detect_stat_seq_kcusum = zeros(sample_size,1);
    


    for t_idx = 2:2:sample_size


        post_change_block = post_change_sample((t_idx-1):t_idx,:);
        pre_change_block = pre_change_sample(randperm(sample_size_ref,2),:);
        
        x0 = pre_change_block(1,:);
        x1 = pre_change_block(2,:);
        y0 = post_change_block(1,:);
        y1 = post_change_block(2,:);
        
        tmp_v = h_RBF(x0,x1,y0,y1,sgma) - delta;


        detect_stat_seq_kcusum(t_idx) = max(0,detect_stat_seq_kcusum(t_idx-1) + tmp_v);

    end
    
    if mod(sample_size,2) == 0
        detect_stat_seq_kcusum(3:2:(sample_size-1)) = detect_stat_seq_kcusum(2:2:(sample_size-2));
    else
        detect_stat_seq_kcusum(3:2:(sample_size-2)) = detect_stat_seq_kcusum(2:2:(sample_size-3));
    end
    
%     detect_stat_seq_cusum = zeros(sample_size,1);
%     
%     for t_idx = 1:sample_size
% 
% 
%         x0 = post_change_sample(t_idx,:);
% 
%         tmp_v = 
%         
%         
%         MMD_u_sq(x0,x1,sgma) + MMD_u_sq(y0,y1,sgma) - MMD_u_sq(x0,y1,sgma) - MMD_u_sq(x1,y0,sgma) - delta;
% 
% 
%         detect_stat_seq_kcusum(t_idx) = max(0,detect_stat_seq_kcusum(t_idx-1) + tmp_v);
% 
%     end    
    
    detect_stat_seq = detect_stat_seq_kcusum;
    
          
end