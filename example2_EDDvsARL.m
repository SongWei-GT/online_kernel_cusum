clear all
% This code can generate the fist panel in Figure 4 in the paper
rng(2023);
%% generate and save detection statistics under H0 and H1
% this step is time consuming and we save the generated data

% for myseed = 1:40
%     generate_stat_H0(myseed-1)
% end

unit_len = 20;

% generate and save detection statistics under H1
% generate_stat_H1()
 
%% settings
sample_dim = 20;

% pre- and post-change distribution parameters
mean1 = 0;
std1 = sqrt(1);

% post-change is Gaussian mixture
post_change_type = "mixG";
mix_p = 0.3;
mean2 = 0;
std2 = 2;   


filename = [];
setting =  ')std_Gaussian_dim_' +  string(sample_dim) ;
filename = [ filename  setting + 'det_stat_max_H0.csv'];
filename = [ filename  setting + 'det_stat_sliding_window_H0.csv'];
filename = [ filename  setting + 'det_stat_H_T2_H0.csv'];
filename = [ filename  setting + 'det_stat_kcusum_H0.csv'];

filename_collection_H0 = filename;

num_file = length(filename_collection_H0);

filename = [];

setting = 'results/(G_to_'+post_change_type+')mix_p_' +  string(mix_p) + '_mean1_' + string(mean1) + '_mean2_' +  string(mean2) + '_std1_' +  string(std1) + '_std2_' +  string(std2);
filename = [ filename  setting + 'det_stat_max_H1.csv'];
filename = [ filename  setting + 'det_stat_sliding_window_H1.csv'];
filename = [ filename  setting + 'det_stat_H_T2_H1.csv'];
filename = [ filename  setting + 'det_stat_kcusum_H1.csv'];


filename_collection_H1 = filename;

EDD_collection = [];

target_ARL = 10.^([2.5:0.5:5]);
sample_size_true = 2000;
lower_quantile = exp(-sample_size_true./target_ARL);

%% load the data and calculate threshold b and EDD
try
    for file_idx = 1:num_file
        det_stat_max_H0 = zeros(1,unit_len*40);

        sample_dim = 20;
        gamma = 1;
        myidx = 1;

        for tmp_idx = 1:40


            filename = 'results/(seed_' + string(myidx-1) + filename_collection_H0(file_idx);
            det_stat_max_H0((myidx-1)*unit_len+1:(myidx)*unit_len)= readmatrix(filename);


            myidx = myidx + 1;
        end

        max_thres_b = quantile(det_stat_max_H0,lower_quantile);


        b_num = length(max_thres_b);
        rep_num = 800;


        max_DD = zeros(b_num,1);



        filename = filename_collection_H1(file_idx);

        det_stat_max_H1 = readmatrix(filename);


        for b_idx = 1:b_num
            for i = 1:rep_num  


                max_temp = find(det_stat_max_H1(i,:) > max_thres_b(b_idx));


                if length(max_temp) == 0
                    max_DD(b_idx,i) =  inf;
                else
                    max_DD(b_idx,i) =  max_temp(1);
                end



            end
        end


        max_EDD = mean(max_DD,2);
        EDD_collection = [EDD_collection
            max_EDD'];
    end



catch
    'mean2: ' +  string(mean2) + ' std2: ' +  string(std2)
    
end


%% visualization
try
    figure();


    p(1) = plot((log10(target_ARL)),EDD_collection(1,:), '-o','MarkerSize',10,'linewidth',3); hold on;
    p(2) = plot((log10(target_ARL)),EDD_collection(2,:), '-.o','MarkerSize',10,'linewidth',3);  
    p(3) = plot((log10(target_ARL)),EDD_collection(4,:), '-.d','MarkerSize',10,'linewidth',3);  
    p(4) = plot((log10(target_ARL)),EDD_collection(3,:), '-.*','MarkerSize',10,'linewidth',3);  

    set(p(1), 'Color', [204,0,0] ./ 255); 
    set(p(2), 'Color', [255, 153, 51] ./ 255); 
    set(p(3), 'Color', [0, 76, 153] ./ 255); 
    set(p(4), 'Color', [153, 76, 0] ./ 255); 


    MATLAB_title = ('\mu = ' + string(mean2) + ', \sigma^2 = ' + string(std2^2));
    title(MATLAB_title, 'Interpreter', 'Tex','fontsize',18)
    set(gca,'fontsize',18,'LooseInset',get(gca,'TightInset')); grid on
    width = 280*1.5;
    set(gcf,'unit','points','PaperUnits','points','PaperPosition',[0,0,width,width/5*4],...
    'color','w','PaperSize',[width, width/5*4]);
    xlabel('\textbf{log ARL}','fontsize',18,'Interpreter','latex');
    ylabel('\textbf{EDD}','fontsize',18,'Interpreter','latex')
    grid(gca,'minor')
    ylim([0 22])
    xlim([2.95 5.05])
    grid off;

    print('plots/example2.pdf','-dpdf','-fillpage')


catch
    'mean2: ' +  string(mean2) + ' std2: ' +  string(std2)
    
end  


