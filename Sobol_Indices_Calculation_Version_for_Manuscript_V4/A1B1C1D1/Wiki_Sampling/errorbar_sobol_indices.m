% errorbar_sobol_indices.m
% --------------------------------------------
% For sample variance, the degree of freedom is n-1 instead of n
%
%
% Input dir/
% D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1C1D1\Wiki_Sampling\Table

clear;

%% read source table from directories

input_dir = 'D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1C1D1\Wiki_Sampling\Table\';

%% box plot of sobol indices results from 3 trials
trial_num = 10;
q = 4;    

Input_Number = 11;
Output_Number = 9;
    
single_figure = 'ErrorBar_Single_Sobol.svg';
total_figure = 'ErrorBar_Total_Sobol.svg';

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
inputVariable = {'Atrsc','Ptrsc','KdeA','AdeA','PdeA','KpolyA','ApolyA','PpolyA','Kdgrd','Adgrd','Pdgrd'};
outputQuantity = {'mean L/S','amp. of L/S','phase of L/S','mean L+S','amp. of L+S','phase of L+S','mean L','amp. of L','phase of L'};


% box plot for single sobol indices
figure;
for i=1:length(outputQuantity)
    
    tmp = Single_Sobol(:,i,:);
    tmp = reshape(tmp,[Input_Number,trial_num]);
    tmp = tmp';
    
    
    subplot(3,3,i)
    x=1:1:Input_Number;
    y = mean(tmp);
    err = std(tmp);
    errorbar(x,y,err,'LineWidth',1,'LineStyle','None','Marker','.','MarkerSize',5,'CapSize',4);
    
    
    ax = gca;
    ax.TickLength = [0.03,0.08];
    ax.YLim = [0 1];
    ax.LineWidth = 2;
    ax.FontSize = 10;
    ax.XLim  = [0 12];
    xticks(1:1:11)
    set(gca,'xticklabel',{[]});

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
    x=1:1:Input_Number;
    y = mean(tmp);
    err = std(tmp);
    errorbar(x,y,err,'LineWidth',1,'LineStyle','None','Marker','.','MarkerSize',5,'CapSize',4);
    ax.XLim  = [0 12];
    xticks(1:1:11)

    ax = gca;
    ax.TickLength = [0.03,0.08];
    ax.YLim = [0 1];
    ax.LineWidth = 2;
    ax.FontSize = 10;
    set(gca,'xticklabel',{[]});
    
    title(outputQuantity{i});
end
sgtitle('Total Sobol indices');
saveas(gcf,total_figure);


