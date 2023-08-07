% boxplot_sobol_indices.m
% --------------------------------------------
% For sample variance, the degree of freedom is n-1 instead of n
%
%
% Input dir/
% D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1C1D1\Nested_Sampling\Table

clear;

%% read source table from directories

input_dir = 'D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1C1D1\Nested_Sampling\Table\';

%% box plot of sobol indices results from 3 trials
trial_num = 10;
q = 4;    

Input_Number = 11;
Output_Number = 9;
    
single_figure = 'Boxplot_Single_Sobol.png';
total_figure = 'Boxplot_Total_Sobol.png';

Single_Sobol = zeros(Input_Number,Output_Number,trial_num);
Total_Sobol = zeros(Input_Number,Output_Number,trial_num);

for trial=1:trial_num 

    % single sobol indices
    single_filename = strcat(input_dir,'Trial_',num2str(trial),'_','SingleSobol_GroupNumber25000_GroupSize_4_Method_Latin.xlsx');

    single_T = readtable(single_filename,'ReadRowNames',true);
    single_T = table2array(single_T);
    single_T = 1-(1-single_T)*q/(q-1); % fix sample variance
    
    Single_Sobol(:,:,trial) = single_T;

    % total sobol indices
    total_filename =  strcat(input_dir,'Trial_',num2str(trial),'_','TotalSobol_GroupNumber25000_GroupSize_4_Method_Latin.xlsx');

    total_T = readtable(total_filename,'ReadRowNames',true);
    total_T = table2array(total_T);
    total_T = total_T*q/(q-1); % fix sample variance
    
    Total_Sobol(:,:,trial) = total_T;
end

%% plotting stuff   
inputVariable = {'Atrsc','Ptrsc','KdeA','AdeA','PdeA','KpolyA','ApolyA','PpolyA','Kdgrd','Adgrd','Pdgrd'};
outputQuantity = {'meanls','rals','Pls','meanss','rass','Pss','meanL','raL','PL'};


% box plot for single sobol indices
figure;
for i=1:length(outputQuantity)
    
    tmp = Single_Sobol(:,i,:);
    tmp = reshape(tmp,[Input_Number,trial_num]);
    tmp = tmp';
    
    
    subplot(3,3,i)
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
    
    
    subplot(3,3,i)
    boxplot(tmp,'Labels',inputVariable);
    xtickangle(45)
    title(outputQuantity{i});
end
sgtitle('Total Sobol indices');
saveas(gcf,total_figure);


