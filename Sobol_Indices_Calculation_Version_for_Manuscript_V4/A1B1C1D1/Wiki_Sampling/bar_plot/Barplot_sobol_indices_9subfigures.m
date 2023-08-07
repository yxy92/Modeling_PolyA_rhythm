% barplot_sobol_indices_9subfigures.m
% --------------------------------------------
% For sample variance, the degree of freedom is n-1 instead of n
%
%
% Input dir/
% D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1D1\Wiki_Sampling\Table

clear;

%% read source table from directories

input_dir = 'D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1D1\Wiki_Sampling\Table\';

%% Bar plot of sobol indices results from 3 trials
trial_num = 10;
q = 4;    

Input_Number = 8;
Output_Number = 9;
    
single_figure = 'Barplot_Single_Sobol.png';
total_figure = 'Barplot_Total_Sobol.png';

Single_Sobol = zeros(Input_Number,Output_Number,trial_num);
Total_Sobol = zeros(Input_Number,Output_Number,trial_num);

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
end

%% plotting stuff   
inputVariable = {'Atrsc','Ptrsc','KdeA','AdeA','PdeA','Kdgrd','Adgrd','Pdgrd'};
outputQuantity = {'meanls','rals','Pls','meanss','rass','Pss','meanL','raL','PL'};


% Bar plot for single sobol indices
figure;
for i=1:3
    for j=1:3
        
        tmp = Single_Sobol(:,3*(i-1)+j,:);
        tmp = reshape(tmp,[Input_Number,trial_num]);
        tmp = tmp';


        subplot(3,3,3*(i-1)+j);
        barplot(tmp,'Labels',inputVariable);



        ax = gca;
        ax.TickLength = [0.03,0.08];
        ax.YLim = [0 1];
        ax.LineWidth = 0.9;
        ax.FontSize = 10;
        ax.XTickLabelRotation = 45;
        title(outputQuantity{i});
    end
end

saveas(gcf,single_figure);


%% Bar plot for total sobol indices
figure;
for i=1:length(outputQuantity)
    
    tmp = Total_Sobol(:,i,:);
    tmp = reshape(tmp,[Input_Number,trial_num]);
    tmp = tmp';
    
    
    subplot(3,3,i)
    barplot(tmp,'Labels',inputVariable);
    
    ax = gca;
    ax.TickLength = [0.03,0.08];
    ax.YLim = [0 1];
    ax.LineWidth = 0.9;
    ax.FontSize = 10;
    ax.XTickLabelRotation = 45;
 
    
    
    title(outputQuantity{i});
end
sgtitle('Total Sobol indices');
saveas(gcf,total_figure);

