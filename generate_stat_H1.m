function generate_stat_H1()
%This code generates detection statistics under H1
%The pre-change distribution is specified to be the 20-dim standard Gaussian
%The post-change distribution is
    addpath ./funcs

    rep_num = 800;
    
    % block size for Scan B procedure
    Blk_sz = 50;
    % Bmax for online kernel CUSUM
    Blk_sz_max = Blk_sz;
    % search region for  online kernel CUSUM
    omega_B = 2:2:Blk_sz_max;
    % number of per-change blocks
    Num_blk = 15;
    
    sample_size = 50;
    sample_dim = 20;

    % initialize the array to save the detection statistics
    det_stat_H_T2_H1 = zeros(rep_num,sample_size);

    det_stat_OKCUSUM_H1 = zeros(rep_num,sample_size);

    det_stat_kcusum_H1 = zeros(rep_num,sample_size);

    det_stat_ScB_H1 = zeros(rep_num,sample_size);

 
    % generate/load the history data

    pre_change_sample = readmatrix('raw_pre_change_sample_dim20.csv');
    %%
    % median heuristo select the kernel bandwidth

    dist_mat = EuDist2(pre_change_sample);       
    kernel_bandwidth = median(dist_mat(dist_mat ~= 0));

    delta = 1/50;

    % pre- and post-change distribution parameters
    mean1 = 0;
    std1 = sqrt(1);
    
    % post-change is Gaussian mixture
    mix_p = 0.3;
    mean2 = 0;
    std2 = 2;   
    
    for rep_idx = 1:rep_num
        
        post_change_sample_1 = normrnd(mean1,std1,[sample_size sample_dim]);
        post_change_sample_2 = normrnd(mean2,std2,[sample_size sample_dim]);

        post_change_biornd = binornd(1,mix_p,1,sample_size);
        post_change_biornd = logical(post_change_biornd);

        post_change_sample = post_change_sample_2;
        post_change_sample(post_change_biornd,:) = post_change_sample_1(post_change_biornd,:);



        detect_stat_seq = online_kernel_cusum(pre_change_sample,post_change_sample,omega_B,Num_blk,kernel_bandwidth);
        det_stat_OKCUSUM_H1(rep_idx,:) = (detect_stat_seq(:,1));
        det_stat_ScB_H1(rep_idx,:) = (detect_stat_seq(:,2));



        det_stat_H_T2_H1(rep_idx,:) = (transpose(H_T2_detection_statistic(pre_change_sample,post_change_sample,Blk_sz)));
        det_stat_kcusum_H1(rep_idx,:) = (cusum_detection_statistic(pre_change_sample,post_change_sample,kernel_bandwidth, delta));


    end


    setting = 'results/(G_to_mixG)mix_p_' +  string(mix_p) + '_mean1_' + string(mean1) + '_mean2_' +  string(mean2) + '_std1_' +  string(std1) + '_std2_' +  string(std2);
    filename = setting + 'det_stat_max_H1.csv';
    writematrix(det_stat_OKCUSUM_H1,filename)
    filename = setting + 'det_stat_sliding_window_H1.csv';
    writematrix(det_stat_ScB_H1,filename)
    filename = setting + 'det_stat_H_T2_H1.csv';
    writematrix(det_stat_H_T2_H1,filename)
    filename = setting + 'det_stat_kcusum_H1.csv';
    writematrix(det_stat_kcusum_H1,filename)


end
