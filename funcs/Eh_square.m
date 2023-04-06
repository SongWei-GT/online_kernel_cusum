function [Eh_sq] = Eh_square(pre_change_sample,kernel_bandwidth)
    [sz_ref,dim_ref]=size(pre_change_sample);

    sample_size = floor(sz_ref/4);
    
    x = pre_change_sample(1:sample_size,:);
    x_p = pre_change_sample((sample_size+1):2*sample_size,:);

    y = pre_change_sample((2*sample_size+1):3*sample_size,:);
    y_p = pre_change_sample((3*sample_size+1):4*sample_size,:);


    D1 = diag(EuDist2(x,x_p,0));
    
    D2 = diag(EuDist2(y,y_p,0));
    
    D3 = diag(EuDist2(x,y_p,0));
    
    D4 = diag(EuDist2(x_p,y,0));
    
    
    K = (exp(-D1/(2*kernel_bandwidth.^2)) + exp(-D2/(2*kernel_bandwidth.^2)) - exp(-D3/(2*kernel_bandwidth.^2)) - exp(-D4/(2*kernel_bandwidth.^2))).^2;

    Eh_sq = mean(K);             
end