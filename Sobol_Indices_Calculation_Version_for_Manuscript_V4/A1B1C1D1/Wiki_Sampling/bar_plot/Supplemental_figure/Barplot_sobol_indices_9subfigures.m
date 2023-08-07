% barplot_sobol_indices_9subfigures.m
% -----------------------------------------
% Input dir/
% D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1C1D1\Wiki_Sampling\Table

clear;

%% read source table from directories

input_dir = 'D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1C1D1\Wiki_Sampling\Table\';

%% Bar plot of sobol indices results from 3 trials
trial_num = 10;
q = 4;    

Input_Number = 11;
Output_Number = 9;
   
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
inputVariable = {'Atrsc','Ptrsc','Kdgrd','Adgrd','Pdgrd','KdeA','AdeA','PdeA','KpolyA','ApolyA','PpolyA'};
outputQuantity = {'meanls','rals','Pls','meanss','rass','Pss','meanL','raL','PL'};


% calculate mean value and s.d for bar plot
Single_Sobol_mean = mean(Single_Sobol,3);
Total_Sobol_mean =  mean(Total_Sobol,3);

Single_Sobol_sd =  std(Single_Sobol,[],3);
Total_Sobol_sd =  std(Total_Sobol,[],3);

% flip position of Dgrd and DeA process such that the input variables are
% Trsc, Dgrd, DeA, PolyA
Single_Sobol_mean = Single_Sobol_mean([1 2 9 10 11 3 4 5 6 7 8],:);
Total_Sobol_mean = Total_Sobol_mean([1 2 9 10 11 3 4 5 6 7 8],:);

Single_Sobol_sd = Single_Sobol_sd([1 2 9 10 11 3 4 5 6 7 8],:);
Total_Sobol_sd = Total_Sobol_sd([1 2 9 10 11 3 4 5 6 7 8],:);




% color scheme
color1 = [73 0 146]/255;  % Trsc
color2 = [0 109 219]/255; % Dgrd
color3 = [146 73 0]/255;  % DeA
color4 = [219 109 0]/255;     % polyA


% for single sobol
figure;
for i=1:Output_Number
    % plot for single sobol, ylim= [0 0.5]
    subplot(3,3,i);
    % bar graph
    x=1:Input_Number;
    hB =bar(x,Single_Sobol_mean(:,i)');
    hB.FaceColor = 'flat';
    for k=1:2
        hB.CData(k,:)= color1;
    end
    for k=3:5
        hB.CData(k,:) = color2;
    end
    for k=6:8
        hB.CData(k,:) = color3;
    end
    for k=9:11
        hB.CData(k,:) = color4;
    end

    % customize axes
    ax = gca;
    ax.TickLength = [0.02,0.02];
    ax.YLim = [0 1];
    ax.YTick = [0 0.5 1];
    ax.LineWidth = 1;
    ax.FontSize = 12;    
end
sgtitle('Single Sobol');
saveas(gcf,'Barplot_Single_Sobol_9_subplots.png');

% plot for total sobol ,ylim = [0 1] 
figure;
for i=1:Output_Number  
    subplot(3,3,i);
    % bar graph
    x=1:Input_Number;
    hB =bar(x,Total_Sobol_mean(:,i)');
    hB.FaceColor = 'flat';
    for k=1:2
        hB.CData(k,:)= color1;
    end
    for k=3:5
        hB.CData(k,:) = color2;
    end
    for k=6:8
        hB.CData(k,:) = color3;
    end
    for k=9:11
        hB.CData(k,:) = color4;
    end
    % customize axes
    ax = gca;
    ax.TickLength = [0.02,0.02];
    ax.YLim = [0 1];
    ax.YTick = [0 0.5 1];
    ax.LineWidth = 1;
    ax.FontSize = 12;
    
  
end 
sgtitle('Total Sobol');
saveas(gcf,'Barplot_Total_Sobol_9_subplots.png');




