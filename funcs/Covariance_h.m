function [Covh] = Covariance_h(pre_change_sample,kernel_bandwidth)
    [sz_ref,dim_ref]=size(pre_change_sample);

    sample_size = floor(sz_ref/6);
    
    x = pre_change_sample(1:sample_size,:);
    x_p = pre_change_sample((sample_size+1):2*sample_size,:);
    x_pp = pre_change_sample((4*sample_size+1):5*sample_size,:);
    x_ppp = pre_change_sample((5*sample_size+1):6*sample_size,:);

    y = pre_change_sample((2*sample_size+1):3*sample_size,:);
    y_p = pre_change_sample((3*sample_size+1):4*sample_size,:);


    D1 = diag(EuDist2(x,x_p,0));
    
    D2 = diag(EuDist2(y,y_p,0));
    
    D3 = diag(EuDist2(x,y_p,0));
    
    D4 = diag(EuDist2(x_p,y,0));
    
    D1_p = diag(EuDist2(x_pp,x_ppp,0));
    
    D2_p = diag(EuDist2(y,y_p,0));
    
    D3_p = diag(EuDist2(x_pp,y_p,0));
    
    D4_p = diag(EuDist2(x_ppp,y,0));
    
    
    K1 = (exp(-D1/(2*kernel_bandwidth.^2)) + exp(-D2/(2*kernel_bandwidth.^2)) - exp(-D3/(2*kernel_bandwidth.^2)) - exp(-D4/(2*kernel_bandwidth.^2)));
    K2 = (exp(-D1_p/(2*kernel_bandwidth.^2)) + exp(-D2_p/(2*kernel_bandwidth.^2)) - exp(-D3_p/(2*kernel_bandwidth.^2)) - exp(-D4_p/(2*kernel_bandwidth.^2)));
    

    
    Covh = mean(K1.*K2) - mean(K1) * mean(K2);

end