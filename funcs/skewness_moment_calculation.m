sample_size_true = 10000;
% sample_dim = 20;
% 
% x_0 = normrnd(0,1,[sample_size_true sample_dim]);
% x_1 = normrnd(0,1,[sample_size_true sample_dim]);
% x_2 = normrnd(0,1,[sample_size_true sample_dim]);
% x_3 = normrnd(0,1,[sample_size_true sample_dim]);
% x_4 = normrnd(0,1,[sample_size_true sample_dim]);
% x_5 = normrnd(0,1,[sample_size_true sample_dim]);
% y_0 = normrnd(0,1,[sample_size_true sample_dim]);
% y_1 = normrnd(0,1,[sample_size_true sample_dim]);
% y_2 = normrnd(0,1,[sample_size_true sample_dim]);


sample_dim = 1;
x_0 = laprnd(sample_size_true, sample_dim, 0, 1);
x_1 = laprnd(sample_size_true, sample_dim, 0, 1);
x_2 = laprnd(sample_size_true, sample_dim, 0, 1);
x_3 = laprnd(sample_size_true, sample_dim, 0, 1);
x_4 = laprnd(sample_size_true, sample_dim, 0, 1);
x_5 = laprnd(sample_size_true, sample_dim, 0, 1);
y_0 = laprnd(sample_size_true, sample_dim, 0, 1);
y_1 = laprnd(sample_size_true, sample_dim, 0, 1);
y_2 = laprnd(sample_size_true, sample_dim, 0, 1);


% x_0 = exprnd(1,[sample_size_true sample_dim]);
% x_1 = exprnd(1,[sample_size_true sample_dim]);
% x_2 = exprnd(1,[sample_size_true sample_dim]);
% x_3 = exprnd(1,[sample_size_true sample_dim]);
% x_4 = exprnd(1,[sample_size_true sample_dim]);
% x_5 = exprnd(1,[sample_size_true sample_dim]);
% y_0 = exprnd(1,[sample_size_true sample_dim]);
% y_1 = exprnd(1,[sample_size_true sample_dim]);
% y_2 = exprnd(1,[sample_size_true sample_dim]);

%variance calculation
Eh1h2h3_1 = 0;
Eh1h2h3_2 = 0;
Eh1h2h3_3 = 0;
Eh1h2h3_4 = 0;
Eh_cube = 0;
Eh1sqh2 = 0;
Eh1h2 = 0;
Eh_sq = 0;
Eh1 = 0;
Eh2 = 0;
%Eh = zeros(1,sample_size);

% centered_x = x_0(1:10000,:) - mean(x_0(1:10000,:),1);
% kernel_bandwidth = median(sqrt(sum(centered_x.^2,2)));

dist_mat = EuDist2(x_0);       
kernel_bandwidth = median(dist_mat(dist_mat ~= 0));

Eh = zeros(1,sample_size);

for i = 1:sample_size_true
    temp_0_1_0_1 = h_RBF(x_0(i,:),x_1(i,:),y_0(i,:),y_1(i,:),kernel_bandwidth);
    temp_1_2_1_2 = h_RBF(x_1(i,:),x_2(i,:),y_1(i,:),y_2(i,:),kernel_bandwidth);
    temp_2_0_2_0 = h_RBF(x_2(i,:),x_0(i,:),y_2(i,:),y_0(i,:),kernel_bandwidth);
    temp_3_4_2_0 = h_RBF(x_3(i,:),x_4(i,:),y_2(i,:),y_0(i,:),kernel_bandwidth);
    temp_2_3_1_2 = h_RBF(x_2(i,:),x_3(i,:),y_1(i,:),y_2(i,:),kernel_bandwidth);
    temp_4_5_2_0 = h_RBF(x_4(i,:),x_5(i,:),y_2(i,:),y_0(i,:),kernel_bandwidth);
    temp_2_3_0_1 = h_RBF(x_2(i,:),x_3(i,:),y_0(i,:),y_1(i,:),kernel_bandwidth);
    temp_4_5_0_1 = h_RBF(x_4(i,:),x_5(i,:),y_0(i,:),y_1(i,:),kernel_bandwidth);
    Eh_sq = Eh_sq + temp_0_1_0_1.^2./sample_size_true; 
    Eh1h2 = Eh1h2 + temp_2_3_0_1.*temp_0_1_0_1./sample_size_true; 
    Eh1 = Eh1 + temp_0_1_0_1./sample_size_true; 
    Eh2 = Eh2 + temp_2_3_0_1./sample_size_true; 
    Eh1h2h3_1 = Eh1h2h3_1 + (temp_0_1_0_1.*temp_1_2_1_2.*temp_2_0_2_0)./sample_size_true; 
    Eh1h2h3_2 = Eh1h2h3_2 + (temp_0_1_0_1.*temp_1_2_1_2.*temp_3_4_2_0)./sample_size_true; 
    Eh1h2h3_3 = Eh1h2h3_3 + (temp_0_1_0_1.*temp_2_3_1_2.*temp_4_5_2_0)./sample_size_true; 
    Eh1h2h3_4 = Eh1h2h3_4 + (temp_0_1_0_1.*temp_2_3_0_1.*temp_4_5_0_1)./sample_size_true; 
    Eh_cube = Eh_cube + temp_0_1_0_1.^3./sample_size_true; 
    Eh1sqh2 = Eh1sqh2 + temp_0_1_0_1.^2.*temp_2_3_0_1./sample_size_true; 
end  

Covh = Eh1h2 - Eh1.*Eh2;

variance_est = zeros(1,length(omega_B));


for i = 1:length(omega_B)
    variance_est(i) = (Eh_sq./Num_blk + (1-1./Num_blk).*Covh)./nchoosek(omega_B(i),2);   
end


%c3B_lap = zeros(1,length(omega_B));
c3B_lap = zeros(1,length(omega_B));

B_idx = 0;
for B = omega_B
    B_idx = B_idx + 1;
    temp_term_1 = (Eh1h2h3_1 + Eh1h2h3_2.*3.*(Num_blk-1) + Eh1h2h3_3.*(Num_blk-1).*(Num_blk-2))./Num_blk.^2;
    temp_term_2 = (Eh_cube + Eh1sqh2.*3.*(Num_blk-1) + Eh1h2h3_4.*(Num_blk-1).*(Num_blk-2))./Num_blk.^2;
    temp_term = 8.*(B-2)./B.^2./(B-1).^2.*temp_term_1 + 4./B.^2./(B-1).^2.*temp_term_2;
    c3B_lap(B_idx) = temp_term./(variance_est(B_idx).^(3./2));
end