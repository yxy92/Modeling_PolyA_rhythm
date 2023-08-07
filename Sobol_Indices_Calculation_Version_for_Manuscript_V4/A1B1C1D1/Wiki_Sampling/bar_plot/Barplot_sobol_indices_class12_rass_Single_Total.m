% Barplot_sobol_indices_class12_rass_Single_Total.m
% --------------------------------------------
% For sample variance, the degree of freedom is n-1 instead of n
%
%
% Input dir/
% D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1D1\Wiki_Sampling\Table

clear;

%% read source table from directories

input_dir = 'D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1C1D1\Wiki_Sampling\Table\';
output_dir = 'D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1C1D1\Wiki_Sampling\bar_plot\Main_Text_Figure\';


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
inputVariable = {'Atrsc','Ptrsc','KdeA','AdeA','PdeA','Kdgrd','Adgrd','Pdgrd'};
outputQuantity = {'meanls','rals','Pls','meanss','rass','Pss','meanL','raL','PL'};


% color scheme
color1 = [73 0 146]/255;  % Trsc
color2 = [0 109 219]/255; % Dgrd
color3 = [146 73 0]/255;  % DeA
color4 = [219 109 0]/255; % polyA


% calculate mean value and s.d for bar plot
Single_Sobol_mean = mean(Single_Sobol,3);
Total_Sobol_mean =  mean(Total_Sobol,3);

Single_Sobol_sd =  std(Single_Sobol,[],3);
Total_Sobol_sd =  std(Total_Sobol,[],3);

% flip position of Dgrd and DeA process such that the input variables are
% Trsc, Dgrd, DeA, polyA
Single_Sobol_mean = Single_Sobol_mean([1 2 9 10 11 3 4 5 6 7 8],:);
Total_Sobol_mean = Total_Sobol_mean([1 2 9 10 11 3 4 5 6 7 8],:);



Single_Sobol_sd = Single_Sobol_sd([1 2 9 10 11 3 4 5 6 7 8],:);
Total_Sobol_sd = Total_Sobol_sd([1 2 9 10 11 3 4 5 6 7 8],:);


i= 5; %rass
% plot for single sobol, ylim= [0 0.5]
figure;
pos = get(gca, 'Position');
pos(1) = 0.2;
pos(2) = 0.3;
pos(3) = 0.7;
pos(4) = 0.5;
set(gca, 'Position', pos)

% bar graph
%x=[1 2 4 5 6 8 9 10 12 13 14];   % position of each bar
% mix single with total sobol for grouped bar plot
single_tmp = Single_Sobol_mean(:,i);
total_tmp = Total_Sobol_mean(:,i);
Sobol_mean = cat(2,single_tmp,total_tmp);


hB = bar(Sobol_mean,'BarWidth',1,'FaceColor','flat');


% customize face color of each bar
for k=1:2
    for p=1:2
        hB(k).CData(p,:)= color1;
    end
    for p=3:5
        hB(k).CData(p,:) = color2;
    end
    for p=6:8
        hB(k).CData(p,:) = color3;
    end
    for p=9:11
        hB(k).CData(p,:) = color4;
    end   
end

% cumstomize text on top of each bar
for j=1:Input_Number
    text(j-0.2,Sobol_mean(j,1),'S','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20); % Single
    text(j+0.2,Sobol_mean(j,2),'T','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20); % Total
end
% customize axes
ax = gca;
ax.TickLength = [0.01,0.01];
ax.YLim = [0 0.5];
%ax.XTick = [0.5 3 4.5 7];
ax.YTick = [0 0.5 1];
ax.LineWidth = 2;
ax.FontSize = 40;

saveas(gcf,strcat(output_dir,'Single_Total_Sobol_Group_',outputQuantity{i},'.svg'));
    
    
 




