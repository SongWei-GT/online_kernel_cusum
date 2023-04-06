function [stat_H2] = H_T2_detection_statistic(data_ref,data_test,m0)

    data = [data_ref;data_test];
    T_test = length(data_test);
    T_ref = length(data_ref);


    hatmu0 = mean(data_ref);   hatsigma0 = cov(data_ref);

    stat_H2 = zeros(T_test,1);
    for t = 1:T_test
        jj = T_ref+t;
        stat_H2(t) = (mean(data(jj-m0+1:jj,:))-hatmu0)*inv(hatsigma0)*(mean(data(jj-m0+1:jj,:))-hatmu0)';
    end
end