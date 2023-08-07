% boxplot_sobol_indices.m
% --------------------------------------------
% For sample variance, the degree of freedom is n-1 instead of n
%
%
% Input dir/
% D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1D1\Wiki_Sampling\Table

clear;

%% read source table from directories

input_dir = 'D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\Sarah_Model\Wiki_Sampling\Table\';

%% box plot of sobol indices results from 3 trials
trial_num = 10;
    

Input_Number = 5;
Output_Number = 3;
    
single_figure = 'Boxplot_Single_Sobol.png';
total_figure = 'Boxplot_Total_Sobol.png';

Single_Sobol = zeros(Input_Number,Output_Number,trial_num);
Total_Sobol = zeros(Input_Number,Output_Number,trial_num);

for trial=1:trial_num 

    % single sobol indices
    single_filename = strcat(input_dir,'Trial_',num2str(trial),'SingleSobol_SampleSize_10000Method_Latin.xlsx');

    single_T = readtable(single_filename,'ReadRowNames',true);
    single_T = table2array(single_T);
    
    Single_Sobol(:,:,trial) = single_T;

    % total sobol indices
    total_filename =  strcat(input_dir,'Trial_',num2str(trial),'TotalSobol_SampleSize_10000Method_Latin.xlsx');

    total_T = readtable(total_filename,'ReadRowNames',true);
    total_T = table2array(total_T);
    
    Total_Sobol(:,:,trial) = total_T;
end

%% plotting stuff   
inputVariable = {'Atrsc','Ptrsc','Kdgrd','Adgrd','Pdgrd'};
outputQuantity = {'meanss','rass','Pss'};


% box plot for single sobol indices
figure;
for i=1:length(outputQuantity)
    
    tmp = Single_Sobol(:,i,:);
    tmp = reshape(tmp,[Input_Number,trial_num]);
    tmp = tmp';
    
    
    subplot(3,1,i)
    boxplot(tmp,'Labels',inputVariable);
    xtickangle(45)
    title(outputQuantity{i});
end
sgtitle('Single Sobol indices');
saveas(gcf,single_figure);


% box plot for total sobol indices
figure;
for i=1:length(outputQuantity)
    
    tmp = Total_Sobol(:,i,:);
    tmp = reshape(tmp,[Input_Number,trial_num]);
    tmp = tmp';
    
    
    subplot(3,1,i)
    boxplot(tmp,'Labels',inputVariable);
    xtickangle(45)
    title(outputQuantity{i});
end
sgtitle('Total Sobol indices');
saveas(gcf,total_figure);


