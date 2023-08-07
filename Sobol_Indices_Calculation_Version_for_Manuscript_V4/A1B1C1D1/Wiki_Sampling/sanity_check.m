% sanity_check.m
%
% check if Total sobol indices are larger than single sobol indices
%

clear;

input_dir = 'D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1C1D1\Wiki_Sampling\Table\';

%% box plot of sobol indices results from 3 trials
trial_num = 10;
q = 4;    

Input_Number = 11;
Output_Number = 9;
    

Single_Sobol = zeros(Input_Number,Output_Number,trial_num);
Total_Sobol = zeros(Input_Number,Output_Number,trial_num);

sanity = zeros(Input_Number,Output_Number,trial_num);

difference = zeros(Input_Number,Output_Number,trial_num);
for trial=1:trial_num 

    % single sobol indices
    single_filename = strcat(input_dir,'Trial_',num2str(trial),'SingleSobol_SampleSize_100000Method_Latin.xlsx');

    single_T = readtable(single_filename,'ReadRowNames',true);
    single_T = table2array(single_T);
    
    Single_Sobol(:,:,trial) = single_T;

    % total sobol indices
    total_filename =  strcat(input_dir,'Trial_',num2str(trial),'TotalSobol_SampleSize_100000Method_Latin.xlsx');

    total_T = readtable(total_filename,'ReadRowNames',true);
    total_T = table2array(total_T);
    
    Total_Sobol(:,:,trial) = total_T;
    
    
    % compare single sobol indices with total sobol indices
    sanity(:,:,trial) = Total_Sobol(:,:,trial) > Single_Sobol(:,:,trial);
    
    difference(:,:,trial) = (Total_Sobol(:,:,trial)- Single_Sobol(:,:,trial))./(Total_Sobol(:,:,trial)+Single_Sobol(:,:,trial));
end


%% plotting stuff
inputVariable = {'Atrsc','Ptrsc','KdeA','AdeA','PdeA','KpolyA','ApolyA','PpolyA','Kdgrd','Adgrd','Pdgrd'};
outputQuantity = {'meanls','rals','Pls','meanss','rass','Pss','meanL','raL','PL'};


% box plot for single sobol indices
figure;
for i=1:3
    for j=1:3
    
        tmp = difference(:,3*(i-1)+j,:);
        tmp = reshape(tmp,[Input_Number,trial_num]);
        tmp = tmp';


        subplot(3,3,3*(i-1)+j)
        boxplot(tmp,'Labels',inputVariable);
        pos = get(gca, 'Position');
        pos(1) = 0.08+0.3*(j-1);
        pos(2) = 1.05-0.31*i;
        pos(3) = 0.25;
        pos(4) = 0.135;
        set(gca, 'Position', pos)

        ax = gca;
        ax.TickLength = [0.03,0.08];
        ax.YLim = [0 1];
        ax.LineWidth = 0.9;
        ax.FontSize = 9;
        ax.XTickLabelRotation = 60;        

        title(outputQuantity{3*(i-1)+j});
    end
end
sgtitle('Difference between Total and Single Sobol indices');
saveas(gcf,'boxplot_Difference_Total_Single.png');

