function [detect_stat_seq] = online_kernel_cusum(pre_change_sample,post_change_sample,omega_B,Num_blk,kernel_bandwidth)

    %Detection statistics calculation for online kernel CUSUM
    %This function also outputs the Scan B statistics as a by-product
    %
    %   pre_change_sample: reference/history data, a M-by-d matrix
    %   
    %   post_change_sample： sequential data, a T-by-d matrix
    %   
    %   omega_B: the set for candidate block sizes 
    %   (typically {2,...,Bmax})
    %   
    %   Num_blk: the number of pre-change blocks
    %   (typically we have Num_blk * max(omega_B) < M)
    %
    %   In this implementation we adopt Gaussian RBF kernel
    %   kernel_bandwidth: the kernel bandwidth parameter
    %   
    %   version 1.0 --November/2005
    %
    %   Written by Song Wei (song.wei AT gatech.edu)

    

    
    Eh_sq = Eh_square(pre_change_sample,kernel_bandwidth);
    Covh = Covariance_h(pre_change_sample,kernel_bandwidth);

    variance_est = zeros(1,length(omega_B));


    for i = 1:length(omega_B)
        variance_est(i) = (Eh_sq./Num_blk + (1-1./Num_blk).*Covh)./nchoosek(omega_B(i),2);   
    end
    
    
    B_max = max(omega_B);


    all_sample = [pre_change_sample
    post_change_sample];


    [sample_size_ref,~] = size(pre_change_sample);
    [sample_size,~] = size(post_change_sample);

    reference_data_1 = all_sample((sample_size_ref-(Num_blk+1)*B_max + 2):(sample_size_ref-B_max+1),:);

    time_horizon = (sample_size_ref+1):(sample_size_ref + sample_size);

    detect_stat_seq_max = zeros(sample_size,1);
    detect_stat_seq_sw = zeros(sample_size,1);

    Dxx = zeros(B_max,B_max);
    Dyy_collection = zeros(Num_blk,B_max,B_max);
    Dxy_collection = zeros(Num_blk,B_max,B_max);

    for blk_idx = 1:Num_blk
        pre_change_data = reference_data_1(((blk_idx-1)*B_max + 1):(blk_idx*B_max ),:);
        tmp_Dyy = EuDist2(pre_change_data,[],0);
        Dyy_collection(blk_idx,:,:) = tmp_Dyy;
    end

    t_idx = 0;
    for t = time_horizon

        t_idx = t_idx + 1;


        post_change_block = all_sample((t-B_max+1):t,:);



        if t_idx == 1

            Dxx = EuDist2(post_change_block,post_change_block,0);

            for blk_idx = 1:Num_blk
                pre_change_data = reference_data_1(((blk_idx-1)*B_max + 1):(blk_idx*B_max ),:);
                tmp_Dxy = EuDist2(post_change_block,pre_change_data,0);
                Dxy_collection(blk_idx,:,:) = tmp_Dxy;
            end

        else

            new_sample = post_change_block(end:end,:);
            past_post_sample = post_change_block(1:(end-1),:);

            re_used_Dxx = Dxx(2:end,2:end);
            new_Dxx_off_diag = EuDist2(new_sample,past_post_sample,0);
            new_Dxx_diag = EuDist2(new_sample,new_sample,0);

            Dxx(1:(end-1),1:(end-1)) = re_used_Dxx;
            Dxx(end:(end),end:(end)) = new_Dxx_diag;
            Dxx(end:end,1:(end-1)) = new_Dxx_off_diag;
            Dxx(1:(end-1),end:end) = new_Dxx_off_diag';

            for blk_idx = 1:Num_blk
                pre_change_data = reference_data_1(((blk_idx-1)*B_max + 1):(blk_idx*B_max ),:);

                re_used_Dxy = Dxy_collection(blk_idx,2:end,1:end);
                new_Dxy = EuDist2(new_sample,pre_change_data,0);

                Dxy_collection(blk_idx,1:(end-1),:) = re_used_Dxy;
                Dxy_collection(blk_idx,end:end,:) = new_Dxy;

            end

        end

        max_kernel_CPD_stat = -inf;


        tmp_block_index = 0;

        for B_sz = omega_B

            tmp_block_index = tmp_block_index + 1;

            temp_max_test_stat = 0;

            %temp1 = exp(-Dxx(1:B_sz,1:B_sz)/(2*kernel_bandwidth.^2));
            temp1 = exp(-Dxx((end - B_sz + 1):end,(end - B_sz + 1):end)/(2*kernel_bandwidth.^2));

            for blk_idx = 1:Num_blk

                Dxy = reshape(Dxy_collection(blk_idx,(end - B_sz + 1):end,(end - B_sz + 1):end),B_sz,B_sz);
                Dyy = reshape(Dyy_collection(blk_idx,(end - B_sz + 1):end,(end - B_sz + 1):end),B_sz,B_sz);             
                temp3 = exp(-Dxy/(2*kernel_bandwidth.^2));
                temp2 = exp(-Dyy/(2*kernel_bandwidth.^2));
                temp = temp1 + temp2 - temp3 - temp3';

                MMD2 = sum(2.*sum(temp - triu(temp)))./B_sz./(B_sz-1);  
                
                
                temp_max_test_stat  = temp_max_test_stat + MMD2./Num_blk;

            end

            temp_max_test_stat = temp_max_test_stat./sqrt(variance_est(tmp_block_index));

            if temp_max_test_stat > max_kernel_CPD_stat
                max_kernel_CPD_stat = temp_max_test_stat;
            end
        end

        detect_stat_seq_max(t_idx) = max_kernel_CPD_stat;
        detect_stat_seq_sw(t_idx) = temp_max_test_stat;
    end    
    
    detect_stat_seq = [detect_stat_seq_max detect_stat_seq_sw];
            
end