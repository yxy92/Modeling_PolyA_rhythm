% Barplot_sobol_indices_3_figure_for_phase.m
% --------------------------------------------
% For sample variance, the degree of freedom is n-1 instead of n
%
%
% Input dir/
% D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1D1\Wiki_Sampling\Table

clear;

%% read source table from directories

input_dir = 'D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1D1\Wiki_Sampling\Table\';
output_dir = 'D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1D1\Wiki_Sampling\bar_plot\Main_Text_Figure\';


%% Bar plot of sobol indices results from 3 trials
trial_num = 10;
q = 4;    

Input_Number = 8;
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
inputVariable = {'Atrsc','Ptrsc','KdeA','AdeA','PdeA','Kdgrd','Adgrd','Pdgrd'};
outputQuantity = {'meanls','rals','Pls','meanss','rass','Pss','meanL','raL','PL'};


% color scheme
color1 = [73 0 146]/255;  % Trsc
color2 = [0 109 219]/255; % Dgrd
color3 = [146 73 0]/255;  % DeA



% calculate mean value and s.d for bar plot
Single_Sobol_mean = mean(Single_Sobol,3);
Total_Sobol_mean =  mean(Total_Sobol,3);

Single_Sobol_sd =  std(Single_Sobol,[],3);
Total_Sobol_sd =  std(Total_Sobol,[],3);

% flip position of Dgrd and DeA process such that the input variables are
% Trsc, Dgrd, DeA
Single_Sobol_mean = Single_Sobol_mean([1 2 6 7 8 3 4 5],:);
Total_Sobol_mean = Total_Sobol_mean([1 2 6 7 8 3 4 5],:);

Single_Sobol_sd = Single_Sobol_sd([1 2 6 7 8 3 4 5],:);
Total_Sobol_sd = Total_Sobol_sd([1 2 6 7 8 3 4 5],:);

for i=3:3:9
    % plot for single sobol, ylim= [0 0.5]
    figure;
    pos = get(gca, 'Position');
    pos(1) = 0.2;
    pos(2) = 0.3;
    pos(3) = 0.7;
    pos(4) = 0.5;
    set(gca, 'Position', pos)
    
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

    % customize axes
    ax = gca;
    ax.TickLength = [0.01,0.01];
    ax.YLim = [0 0.5];
    ax.YTick = [0 0.5 1];
    ax.LineWidth = 2;
    ax.FontSize = 24;
    
    saveas(gcf,strcat(output_dir,'Single_Sobol',outputQuantity{i},'.svg'));
    
    % plot for total sobol ,ylim = [0 1] 
  
    figure;
    pos = get(gca, 'Position');
    pos(1) = 0.2;
    pos(2) = 0.3;
    pos(3) = 0.7;
    pos(4) = 0.5;
    set(gca, 'Position', pos)
    
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

    % customize axes
    ax = gca;
    ax.TickLength = [0.01,0.01];
    ax.YLim = [0 1];
    ax.YTick = [0 0.5 1];
    ax.LineWidth = 2;
    ax.FontSize = 24;
    
    saveas(gcf,strcat(output_dir,'Total_Sobol',outputQuantity{i},'.svg'));
    
 
end




