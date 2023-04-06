clear all
% This code can generate Figure 1 in the paper
% for demonstrating purpose, we only plot the trajectory for a single trial
% here
rng(2023);
%%
addpath ./funcs

% block size for Scan B procedure
Blk_sz = 50;
% Bmax for online kernel CUSUM
Blk_sz_max = Blk_sz;
% search region for  online kernel CUSUM
omega_B = 2:2:Blk_sz_max;
% number of per-change blocks
Num_blk = 15;

% history sample size
sample_size_ref = 10000;
% sequential data sample size
sample_size = 80;

% a very simple setting for illustration
sample_dim = 20;

%%
% pre- and post-change distribution parameters
mean1 = 0;
std1 = sqrt(1);

% post-change is Gaussian mixture
mix_p = 0.3;
mean2 = 0;
std2 = 2;  

% history data generation
pre_change_sample = normrnd(0,1,[sample_size_ref sample_dim]);

% Gaussian RBF kernel bandwidth choice
dist_mat = EuDist2(pre_change_sample);       
kernel_bandwidth = median(dist_mat(dist_mat ~= 0));
%tic

% at the begining, we do not have change in distribution
post_change_sample = normrnd(0,1,[sample_size sample_dim]);

obs_1 = post_change_sample;

% generate detection statistics "under H0"
tic
detect_stat_seq = online_kernel_cusum(pre_change_sample,post_change_sample,omega_B,Num_blk,kernel_bandwidth);
toc


% next consider a distribution shift to Gaussian mixture
post_change_sample_1 = normrnd(mean1,std1,[sample_size sample_dim]);
post_change_sample_2 = normrnd(mean2,std2,[sample_size sample_dim]);

post_change_biornd = binornd(1,mix_p,1,sample_size);
post_change_biornd = logical(post_change_biornd);

post_change_sample = post_change_sample_2;
post_change_sample(post_change_biornd,:) = post_change_sample_1(post_change_biornd,:);


obs_2 = [obs_1
post_change_sample];

% generate detection statistics "under H1"
tic
detect_stat_seq_H1 = online_kernel_cusum(pre_change_sample,post_change_sample,omega_B,Num_blk,kernel_bandwidth);
toc


OKCUSUM_stat = [detect_stat_seq(:,1)
    detect_stat_seq_H1(:,1)];  


ScB_stat = [detect_stat_seq(:,2)
    detect_stat_seq_H1(:,2)];  


figure()

p(1) = plot(OKCUSUM_stat, '-','MarkerSize',10,'linewidth',3.5);  hold on;
p(2) = plot(ScB_stat, '-','MarkerSize',10,'linewidth',3.5);  

set(p(1), 'Color', [204,0,0] ./ 255); 
set(p(2), 'Color', [255, 153, 51] ./ 255); 

set(gca,'fontsize',18,'LooseInset',get(gca,'TightInset')); grid on
width = 350*1.5;
set(gcf,'unit','points','PaperUnits','points','PaperPosition',[0,0,width,width/5],...
'color','w','PaperSize',[width, width/5*3.5]);
xlabel('\textbf{Time t}','fontsize',18,'Interpreter','latex');
ylabel('\textbf{Detection statistic}','fontsize',18,'Interpreter','latex')
hl = legend({'Proposed','Scan $B$'});
set(hl,'Interpreter','latex','FontSize',18,'Location','northwest');
grid(gca,'minor')
grid off;

print('plots/example1.pdf','-dpdf','-fillpage')


