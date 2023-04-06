function generate_stat_H0(myseed)
%This code generates detection statistics under H0
%The pre-change distribution is specified to be the 20-dim standard Gaussian

    addpath ./funcs
    
    rng(myseed+100);
    
    %% settings

    rep_num = 20;
    
    % block size for Scan B procedure
    Blk_sz = 50;
    % Bmax for online kernel CUSUM
    Blk_sz_max = Blk_sz;
    % search region for  online kernel CUSUM
    omega_B = 2:2:Blk_sz_max;
    % number of per-change blocks
    Num_blk = 15;

    % initialize the array to save the detection statistics
    det_stat_H_T2_H0 = zeros(1,rep_num);
    det_stat_OKCUSUM_H0 = zeros(1,rep_num);
    det_stat_kcusum_H0 = zeros(1,rep_num);
    det_stat_ScB_H0 = zeros(1,rep_num);

    % generate/load the history data
    sample_size = 2000;
    sample_dim = 20;
    pre_change_sample = readmatrix('raw_pre_change_sample_dim20.csv');
    
    % median heuristic to select the kernel bandwidth
    dist_mat = EuDist2(pre_change_sample);       
    kernel_bandwidth = median(dist_mat(dist_mat ~= 0));
    
    % KCUSUM hyperparameter
    delta = 1/50;
    
    %% generation

    
    parfor rep_idx = 1:rep_num

        % post-change data also has the same distribution with the history data
        post_change_sample = normrnd(0,1,[sample_size sample_dim]);

        detect_stat_seq = online_kernel_cusum(pre_change_sample,post_change_sample,omega_B,Num_blk,kernel_bandwidth);

        det_stat_OKCUSUM_H0(rep_idx) = max(detect_stat_seq(:,1));
        det_stat_ScB_H0(rep_idx) = max(detect_stat_seq(:,2));

        det_stat_H_T2_H0(rep_idx) = max(transpose(H_T2_detection_statistic(pre_change_sample,post_change_sample,Blk_sz)));


        det_stat_kcusum_H0(rep_idx) = max(cusum_detection_statistic(pre_change_sample,post_change_sample,kernel_bandwidth, delta));


    end
    
    %% save the data

    setting = 'results/(seed_' + string(myseed) + ')std_Gaussian_dim_' +  string(sample_dim) + "_rand_";
    filename = setting + 'det_stat_max_H0.csv';
    writematrix(det_stat_max_H0_rand,filename)
    filename = setting + 'det_stat_sliding_window_H0.csv';
    writematrix(det_stat_sw_rand_H0,filename)



    setting = 'results/(seed_' + string(myseed) + ')std_Gaussian_dim_' +  string(sample_dim) ;
    filename = setting + 'det_stat_max_H0.csv';
    writematrix(det_stat_OKCUSUM_H0,filename)
    filename = setting + 'det_stat_sliding_window_H0.csv';
    writematrix(det_stat_ScB_H0,filename)
    filename = setting + 'det_stat_H_T2_H0.csv';
    writematrix(det_stat_H_T2_H0,filename)
    filename = setting + 'det_stat_kcusum_H0.csv';
    writematrix(det_stat_kcusum_H0,filename)


end