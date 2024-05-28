%%% Ndnf-IN Arousal Tuning in WT and Scn1a+/- mice

% Author: Sophie Liebergall
% Updated: 4/10/23
% Purpose: Analysis of correlation of Ndnf-IN activity with arousal
% measures in mouse

clc;
clear all;
close all;

%% Load in aligned neural and behavioral data

% Set directory where output aligned data files are stored
data_dir1 = 'D:\two_photon\processed\ndnf_cre_arousal';
dir_list=dir([data_dir1 '\*FOV*']); % Create list of files in directory

% Create empty structure to store data
data = struct('name',{},'dF',{},'Fneu',{},'fs',{},'speed',{},'ops',{},...
    'pupil',{},'time',{});

for d_ind=1:length(dir_list)%1:length(dir_list)
    d_ind
    close all;
    %clearvars -except fs cur_dir dir_name dir_list d_ind data_dir save_dir
    
    dir_n=[data_dir1 '\' dir_list(d_ind).name]; % get path to file
    filename = fullfile(dir_n,'Data.mat');
    
    % Load in data file
    if isfile(filename)
        load(filename)
    %else
    %    continue
    end
    
    data(d_ind).name = dir_list(d_ind).name;
    data(d_ind).dF = dF_F0_g;
    data(d_ind).Fneu = Fneu;
    data(d_ind).fs = fs;
    data(d_ind).speed = inst_speed_r;
    data(d_ind).ops = ops;
    data(d_ind).pupil = pupil_dia_f_r;
    data(d_ind).time = tr;
    
end

WT = {'5905','5907','5909','5912'};
scn1a = {'5906','5908','5910','5911','5913'};

% Add mouse column & genotype to data table
for i = 1:length(data)

    % Add mouse column to data table
    mouse = split(data(i).name,'_');
    data(i).mouse = mouse(1);
    
    % Add genotype column to data
    if ismember(mouse(1), WT)
        data(i).genotype = 'WT';
    elseif ismember(mouse(1), scn1a)
        data(i).genotype = 'Scn1a';
    end

end

%% Crop data if pupil recording ends early & remove blinks from pupil data

for i = 1:length(data) % loop over each recording
    
    pupil_data = data(i).pupil;
    % flip pupil data
    reversedPupil = flip(pupil_data);
    % Find the index where the NaN values end
    endIndex = find(~isnan(reversedPupil), 1, 'first');
    % Trim the NaN values at the end
    trimmedPupil = pupil_data(1:end-endIndex+1);
    
    data(i).pupil = trimmedPupil;
    
    if height(data(i).pupil) < width(data(i).dF)
        
        % trim width of other variables
        data(i).dF = data(i).dF(:,1:height(data(i).pupil));
        data(i).speed = data(i).speed(1:height(data(i).pupil));
        data(i).Fneu = data(i).Fneu(:,1:height(data(i).pupil));
        
    end

end

% Remove blinks from pupil data

accel_threshold = 1; % set acceleration threshold for pupil smoothing

for i = 1:length(data) % loop over each recording
    
    data(i).pupil = remove_blinks(data(i).pupil,accel_threshold);

end

%% Calculate skew & kurtosis for each cell in dataset

for i = 1:length(data)
   
    data(i).skew = skewness(data(i).dF, 1, 2);
    data(i).kurt = kurtosis(data(i).dF, 1, 2);
    data(i).range = range(data(i).dF, 2);
    
    deviations = data(i).dF - prctile(data(i).dF,10, 2);
    data(i).rms = sqrt(mean(deviations.^2,2));
    
end

%% Create lists of WT & Scn1a colors for plotting

WT_color = [0/250 0/250 0/250];
scn1a_color = [236/255 28/255 36/255];

%% Calculate cross correlation between neural activity, running data and pupil diameter

for i = 1:length(data) % loop over each recording
    
    clear corr_nr corr_np corr_rp lags_nr lags_np lags_rp...
        corr_NR corr_NP corr_RP lags_NR lags_NP corr_RP zero_corr_NR...
        zero_corr_NP zero_corr_RP corr_npd lags_npd corr_NPD lags_NPD
    
    speed = zscore(data(i).speed);
    pupil = zscore(sgolayfilt(fillmissing(data(i).pupil,'linear'),1,51));
    
    for ii=1:size(data(i).dF,1)
        
        % Data files for each cell for alignment
        neural = zscore(sgolayfilt(data(i).dF(ii,:),1,51));
        
        % calculate cross correlation of neural activity and running
        [corr_nr, lags_nr] = xcorr(neural, speed, 300, 'normalized');
        % calculate cross correlation of neural activity and pupil diameter
        [corr_np, lags_np] = xcorr(neural, pupil, 300, 'normalized');
        % calculate cross correlation of neural activity and pupil diameter
        [corr_npd, lags_npd] = xcorr(neural(1:end-1), diff(pupil), 300, 'normalized');
        
        % Get max correlation and lag for that max correlation value and
        % assign to table
        % Neural activity and running
        [c_NR, l_NR] = max(abs(corr_nr));
        corr_NR(ii) = corr_nr(l_NR); % get max correlation value
        lags_NR(ii) = lags_nr(l_NR); % get lag (index) of max correlation
        zero_c_NR = corr_nr(300); % get zero correlation
        zero_corr_NR(ii) = zero_c_NR; 
        % Neural activity and pupil
        [c_NP, l_NP] = max(abs(corr_np));
        corr_NP(ii) = corr_np(l_NP);
        lags_NP(ii) = lags_np(l_NP);
        zero_c_NP = corr_np(300); % get zero correlation
        zero_corr_NP(ii) = zero_c_NP; 
        % Neural activity and pupil derivative
        [c_NPD, l_NPD] = max(abs(corr_npd));
        corr_NPD(ii) = corr_npd(l_NPD);
        lags_NPD(ii) = lags_npd(l_NPD);
        zero_c_NPD = corr_npd(300); % get zero correlation
        zero_corr_NPD(ii) = zero_c_NPD; 
       
    end
    
    % Get max correlation of running and pupil
    [corr_rp, lags_rp] = xcorr(speed, pupil, 1000, 'normalized');
    [c_RP, l_RP] = max(corr_rp);
    corr_RP = c_RP; % get max correlation value
    lags_RP = lags_rp(l_RP); % get lag (index) of max correlation
    zero_corr_RP = corr_rp(1000); % get zero correlation value
    

    % Add max correlation data for each cell to data table
    data(i).corr_NR = corr_NR;
    data(i).corr_NP = corr_NP;
    data(i).corr_RP = corr_RP;
    data(i).lags_NR = lags_NR;
    data(i).lags_NP = lags_NP;
    data(i).lags_RP = lags_RP;
    data(i).zero_corr_NR = zero_corr_NR;
    data(i).zero_corr_NP = zero_corr_NP;
    data(i).zero_corr_RP = zero_corr_RP;
    data(i).corr_NPD = corr_NPD;
    data(i).lags_NPD = lags_NPD;
    data(i).zero_corr_NPD = zero_corr_NPD;
    
end

%% Calculate shuffled cross correlation between neural activity, running data and pupil diameter

for i = 1:length(data) % loop over each recording
    
    clear corr_nr corr_np lags_nr lags_np ...
        corr_NR corr_NP lags_NR lags_NP zero_corr_NR...
        zero_corr_NP
    
    speed = zscore(flip(data(i).speed));
    pupil = zscore(flip(sgolayfilt(fillmissing(data(i).pupil,'linear'),1,51)));
    
    for ii=1:size(data(i).dF,1)
        
        % Data files for each cell for alignment
        neural = zscore(sgolayfilt(data(i).dF(ii,:),1,51));
        
        % calculate cross correlation of neural activity and running
        [corr_nr, lags_nr] = xcorr(neural, speed, 300, 'normalized');
        % calculate cross correlation of neural activity and pupil diameter
        [corr_np, lags_np] = xcorr(neural, pupil, 300, 'normalized');
        
        % Get max correlation and lag for that max correlation value and
        % assign to table
        % Neural activity and running
        [c_NR, l_NR] = max(abs(corr_nr));
        corr_NR(ii) = corr_nr(l_NR); % get max correlation value
        lags_NR(ii) = lags_nr(l_NR); % get lag (index) of max correlation
        zero_c_NR = corr_nr(300); % get zero correlation
        zero_corr_NR(ii) = zero_c_NR; 
        % Neural activity and pupil
        [c_NP, l_NP] = max(abs(corr_np));
        corr_NP(ii) = corr_np(l_NP);
        lags_NP(ii) = lags_np(l_NP);
        zero_c_NP = corr_np(300); % get zero correlation
        zero_corr_NP(ii) = zero_c_NP; 
       
    end
    
    % Add max correlation data for each cell to data table
    data(i).corr_NR_shuffle = corr_NR;
    data(i).corr_NP_shuffle = corr_NP;
    data(i).lags_NR_shuffle = lags_NR;
    data(i).lags_NP_shuffle = lags_NP;
    data(i).zero_corr_NR_shuffle = zero_corr_NR;
    data(i).zero_corr_NP_shuffle = zero_corr_NP;
    
end

%% Create plots of neural activity, pupil, and speed for each FOV

for i=1:length(data)
    
    % Plot neural activity as colormesh
    subplot(3, 1, 1);
    mesh = pcolor(data(i).dF);
    set(mesh, 'EdgeColor', 'none');
    ylabel('Cell Number');
    
    % Plot running data
    subplot(3, 1, 2);
    plot(data(i).speed, 'k');
    ylabel('Speed');
    xlim([1 width(data(i).dF)]);
    
    % Plot pupil data
    subplot(3, 1, 3);
    plot(data(i).pupil, 'b');
    xlabel('Time');
    ylabel('Pupil Diameter');
    xlim([1 width(data(i).dF)]);
    
    figname = strcat(data(i).name, '_traces.png');
    
    saveas(gca, fullfile(savedir,figname));
    
end

%% Display mean correlation coefficient for each recording

corr_table = table('Size',[0,4],'VariableNames',{'name','NR','NP','RP',},...
    'VariableTypes',{'string','double',...
    'double','double'});

for i = 1:length(data)
    
    corr_table.name(i) = data(i).name;
    corr_table.NR(i) = mean(data(i).corr_NR);
    corr_table.NP(i) = mean(data(i).corr_NP);
    corr_table.RP(i) = mean(data(i).corr_RP);
    corr_table.NPD(i) = mean(data(i).corr_NPD);
    
end


%% Correlation Histograms by Genotype

% Create tables with correlation coefficients for all WT/Scn1a mice

WT_NP = [];
WT_NR = [];
WT_NP_shuffle = [];
WT_NR_shuffle = [];
WT_RP = [];
WT_zero_NP = [];
WT_zero_NR = [];
WT_zero_NP_shuffle = [];
WT_zero_NR_shuffle = [];
WT_zero_RP = [];
WT_NP_lags = [];
WT_NR_lags = [];
WT_NP_lags_shuffle = [];
WT_NR_lags_shuffle = [];
WT_RP_lags = [];
scn1a_NP = [];
scn1a_NR = [];
scn1a_NP_shuffle = [];
scn1a_NR_shuffle = [];
scn1a_RP = [];
scn1a_zero_NP = [];
scn1a_zero_NR = [];
scn1a_zero_NP_shuffle = [];
scn1a_zero_NR_shuffle = [];
scn1a_zero_RP = [];
scn1a_NP_lags = [];
scn1a_NR_lags = [];
scn1a_NP_lags_shuffle = [];
scn1a_NR_lags_shuffle = [];
scn1a_RP_lags = [];

for i=1:length(data) % loop over each recording
        
    if ismember(data(i).mouse(1),WT)
        WT_NP = [WT_NP data(i).corr_NP];
        WT_NR = [WT_NR data(i).corr_NR];
        WT_NP_shuffle = [WT_NP_shuffle data(i).corr_NP_shuffle];
        WT_NR_shuffle = [WT_NR_shuffle data(i).corr_NR_shuffle];
        WT_RP = [WT_RP data(i).corr_RP];
        WT_NP_lags = [WT_NP_lags data(i).lags_NP];
        WT_NR_lags = [WT_NR_lags data(i).lags_NR];
        WT_NP_lags_shuffle = [WT_NP_lags_shuffle data(i).lags_NP_shuffle];
        WT_NR_lags_shuffle = [WT_NR_lags_shuffle data(i).lags_NR_shuffle];
        WT_RP_lags = [WT_RP_lags data(i).lags_RP];
        WT_zero_NP = [WT_zero_NP data(i).zero_corr_NP];
        WT_zero_NR = [WT_zero_NR data(i).zero_corr_NR];
        WT_zero_NP_shuffle = [WT_zero_NP_shuffle data(i).zero_corr_NP_shuffle];
        WT_zero_NR_shuffle = [WT_zero_NR_shuffle data(i).zero_corr_NR_shuffle];
        WT_zero_RP = [WT_zero_RP data(i).zero_corr_RP];
    elseif ismember(data(i).mouse(1),scn1a)
        scn1a_NP = [scn1a_NP data(i).corr_NP];
        scn1a_NR = [scn1a_NR data(i).corr_NR];
        scn1a_NP_shuffle = [scn1a_NP_shuffle data(i).corr_NP_shuffle];
        scn1a_NR_shuffle = [scn1a_NR_shuffle data(i).corr_NR_shuffle];
        scn1a_RP = [scn1a_RP data(i).corr_RP];
        scn1a_NP_lags = [scn1a_NP_lags data(i).lags_NP];
        scn1a_NR_lags = [scn1a_NR_lags data(i).lags_NR];
        scn1a_NP_lags_shuffle = [scn1a_NP_lags_shuffle data(i).lags_NP_shuffle];
        scn1a_NR_lags_shuffle = [scn1a_NR_lags_shuffle data(i).lags_NR_shuffle];
        scn1a_RP_lags = [scn1a_RP_lags data(i).lags_RP];
        scn1a_zero_NP = [scn1a_zero_NP data(i).zero_corr_NP];
        scn1a_zero_NR = [scn1a_zero_NR data(i).zero_corr_NR];
        scn1a_zero_NP_shuffle = [scn1a_zero_NP_shuffle data(i).zero_corr_NP_shuffle];
        scn1a_zero_NR_shuffle = [scn1a_zero_NR_shuffle data(i).zero_corr_NR_shuffle];
        scn1a_zero_RP = [scn1a_zero_RP data(i).zero_corr_RP];
    end
    
end

    
%% Plot Max Correlation Histograms for Each FOV: Neural Activity and Locomotion

%savedir = "C:\Users\liebergals\OneDrive - Children's Hospital of Philadelphia\Goldberg Lab\Sophie\figures\ndnf_scn1a\2P\arousal\running_corr_by_FOV";

for i = 1:width(data)
    
    N = hist(data(i).corr_NR, [-1:0.05:1]); % Create histogram of correlation coeffs
    pd = fitdist(data(i).corr_NR.','normal'); % Fit probability dist to data
    q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
    x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
    y_g = numel(data(i).corr_NR)*0.05*pdf(pd,x_g); % Plot pdf for data
    
    figure;
    b1 = bar([-1:0.05:1],N,'FaceColor','k','EdgeColor','none');
    b1.FaceAlpha = 0.5;
    
    hold on;
    
    plot([mean(data(i).corr_NR) mean(data(i).corr_NR)],[0 max(N)+3],'--','Color','k','LineWidth',2.5);
    hold on;
    txt = 'Mean Corr: ' + string(round(mean(data(i).corr_NR),2));
    text((mean(data(i).corr_NR)-0.05),(max(N)+1),txt,...
        'FontSize',14,'HorizontalAlignment', 'right');
    
    %hold on;plot(x_g,y_g,'g','LineWidth',2);
    %     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
    %     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    t = strcat("Neural Activity and Locomotion ", char(data(i).name));
    title(t,'fontsize',14);
    
    ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
    figname = strcat(data(i).name, '_neural_running_corr.png');
    saveas(gca, fullfile(savedir,figname)); % Save file as png
    
    close;
    
end

%% Plot Max Correlation Histograms for Each FOV: Neural Activity and Pupil

%savedir = "C:\Users\liebergals\OneDrive - Children's Hospital of Philadelphia\Goldberg Lab\Sophie\figures\ndnf_scn1a\2P\arousal\pupil_corr_by_FOV";

for i = 1:width(data)
    
    close all
    
    N = hist(data(i).corr_NP, [-1:0.05:1]); % Create histogram of correlation coeffs
    pd = fitdist(data(i).corr_NP.','normal'); % Fit probability dist to data
    q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
    x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
    y_g = numel(data(i).corr_NP)*0.05*pdf(pd,x_g); % Plot pdf for data
    
    %figure;
    b1 = bar([-1:0.05:1],N,'FaceColor','k','EdgeColor','none');
    b1.FaceAlpha = 0.5;
    
    hold on;
    
    plot([mean(data(i).corr_NP) mean(data(i).corr_NP)],[0 max(N)+3],'--','Color','k','LineWidth',2.5);
    hold on;
    txt = 'Mean Corr: ' + string(round(mean(data(i).corr_NP),2));
    text((mean(data(i).corr_NP)-0.05),(max(N)+1),txt,...
        'FontSize',14,'HorizontalAlignment', 'right');

    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    t = strcat("Neural Activity and Pupil ", char(data(i).name));
    title(t,'fontsize',14);
    
    ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
    figname = strcat(data(i).name, '_neural_pupil_corr.png');
    saveas(gca, fullfile(savedir,figname)); % Save file as png
    
    close;
    
end

%% Plot Max Correlation Histograms of all WT & Scn1a cells: Neural Activity and Locomotion

% WT Neural Activity and Running
N = hist(WT_NR, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(WT_NR.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(WT_NR)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.5;

hold on;

% Scn1a Neural Activity and Running
N = hist(scn1a_NR,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(scn1a_NR.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(scn1a_NR)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.5;

hold on;
plot([mean(WT_NR) mean(WT_NR)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2);
hold on;
plot([mean(scn1a_NR) mean(scn1a_NR)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2);
hold on;
txt_WT = 'WT Mean Corr: ' + string(round(mean(WT_NR),2));
txt_scn1a = 'Scn1a+/- Mean Corr: ' + string(round(mean(scn1a_NR),2));
text((mean(WT_NR)-0.05),(max(N)),txt_WT,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(scn1a_NR)-0.05),(max(N)-5),txt_scn1a,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('Neural Activity and Locomotion','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
% Kolmogorov-Smirnov Test for Normality    
disp('One Sample K-S Testing')
[h, p] = kstest(WT_NR);
disp('WT: ')
disp(h)
disp(p)
[h, p] = kstest(scn1a_NR);
disp('Scn1a+/-: ')
disp(h)
disp(p)

% Kolmogorov-Smirnov Test for Significant Difference Between Distributions

disp('Two Sample K-S Testing')
[h, p] = kstest2(WT_NR, scn1a_NR);
disp(h)
disp(p)
    
%% Plot Max SHUFFLE Correlation Histograms of all WT & Scn1a cells: Neural Activity and Locomotion

% WT Neural Activity and Running
plot_WT_table = WT_NR_shuffle;

N = hist(plot_WT_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_WT_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_WT_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.5;

hold on;

% Scn1a Neural Activity and Running

plot_scn1a_table = scn1a_NR_shuffle;

N = hist(plot_scn1a_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_scn1a_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_scn1a_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.5;

hold on;
plot([mean(plot_WT_table) mean(plot_WT_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(plot_scn1a_table) mean(plot_scn1a_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_WT = 'WT Mean Corr: ' + string(round(mean(plot_WT_table),2));
txt_scn1a = 'Scn1a+/- Mean Corr: ' + string(round(mean(plot_scn1a_table),2));
text((mean(plot_WT_table)-0.05),(max(N)+15),txt_WT,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_scn1a_table)-0.05),(max(N)+5),txt_scn1a,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('Shuffled Neural Activity and Locomotion','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
% Kolmogorov-Smirnov Test for Significant Difference Between Distributions

disp('Two Sample K-S Testing')
[h, p] = kstest2(plot_WT_table, plot_scn1a_table);
disp(h)
disp(p)
    
%% Plot Lag Histograms of all WT & Scn1a cells: Neural Activity and Locomotion

window = 300;

% WT Neural Activity and Running
N = hist(WT_NR_lags, [-window:10:window]); % Create histogram of correlation coeffs
pd = fitdist(WT_NR_lags.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(WT_NR_lags)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-window:10:window],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.5;

hold on;

% Scn1a Neural Activity and Running
N = hist(scn1a_NR_lags,[-window:10:window]); % Create histogram of correlation coeffs
pd = fitdist(scn1a_NR_lags.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(scn1a_NR_lags)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-window:10:window],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.5;

hold on;
plot([mean(WT_NR_lags) mean(WT_NR_lags)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(scn1a_NR_lags) mean(scn1a_NR_lags)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_WT = 'WT Mean Lag: ' + string(round(mean(WT_NR_lags)/fs,2)) + ' s';
txt_scn1a = 'Scn1a+/- Mean Lag: ' + string(round(mean(scn1a_NR_lags)/fs,2)) + ' s';
text((mean(WT_NR)-0.05),(max(N)),txt_WT,'HorizontalAlignment', 'right');
hold on;
text((mean(scn1a_NR)-0.05),(max(N)-5),txt_scn1a,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Lag (s)','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('Lags: Neural Activity and Locomotion','fontsize',14);
    xticks([-window:int16(fs)*2:window]);
    xticklabels([-window:int16(fs)*2:window]/fs);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');

    
%% Plot Shuffled Lag Histograms of all WT & Scn1a cells: Neural Activity and Locomotion

window = 300;

% WT Neural Activity and Running

plot_WT_table = WT_NR_lags_shuffle;

N = hist(plot_WT_table, [-window:10:window]); % Create histogram of correlation coeffs
pd = fitdist(plot_WT_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_WT_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-window:10:window],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.5;

hold on;

% Scn1a Neural Activity and Running

plot_scn1a_table = scn1a_NR_lags_shuffle;

N = hist(plot_scn1a_table,[-window:10:window]); % Create histogram of correlation coeffs
pd = fitdist(plot_scn1a_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_scn1a_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-window:10:window],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.5;

hold on;
plot([mean(plot_WT_table) mean(plot_WT_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(plot_scn1a_table) mean(plot_scn1a_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_WT = 'WT Mean Lag: ' + string(round(mean(plot_WT_table)/fs,2)) + ' s';
txt_scn1a = 'Scn1a+/- Mean Lag: ' + string(round(mean(plot_scn1a_table)/fs,2)) + ' s';
text((mean(plot_WT_table)-0.05),(max(N)),txt_WT,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_scn1a_table)-0.05),(max(N)-5),txt_scn1a,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Lag (s)','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('Shuffled Lags: Neural Activity and Locomotion','fontsize',14);
    xticks([-window:int16(fs)*2:window]);
    xticklabels([-window:int16(fs)*2:window]/fs);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');    

%% Plot Max Correlation Histograms of all WT & Scn1a cells: Neural Activity and Pupil

% WT Neural Activity and Running
N = hist(WT_NP, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(WT_NP.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(WT_NP)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.5;

hold on;

% Scn1a Neural Activity and Running
N = hist(scn1a_NP,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(scn1a_NP.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(scn1a_NP)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.5;

hold on;
plot([mean(WT_NP) mean(WT_NP)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(scn1a_NP) mean(scn1a_NP)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_WT = 'WT Mean Corr: ' + string(round(mean(WT_NP),2));
txt_scn1a = 'Scn1a+/- Mean Corr: ' + string(round(mean(scn1a_NP),2));
text((mean(WT_NP)-0.05),(max(N)),txt_WT,'HorizontalAlignment', 'right');
hold on;
text((mean(scn1a_NP)-0.05),(max(N)-5),txt_scn1a,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('Neural Activity and Pupil','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
% Kolmogorov-Smirnov Test for Significant Difference Between Distributions

disp('Two Sample K-S Testing')
[h, p] = kstest2(WT_NP, scn1a_NP);
disp(h)
disp(p)
    
%% Plot Max SHUFFLE Correlation Histograms of all WT & Scn1a cells: Neural Activity and Pupil

% WT Neural Activity and Pupil
plot_WT_table = WT_NP_shuffle;

N = hist(plot_WT_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_WT_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_WT_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.5;

hold on;

% Scn1a Neural Activity and Pupil

plot_scn1a_table = scn1a_NP_shuffle;

N = hist(plot_scn1a_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_scn1a_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_scn1a_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.5;

hold on;
plot([mean(plot_WT_table) mean(plot_WT_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(plot_scn1a_table) mean(plot_scn1a_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_WT = 'WT Mean Corr: ' + string(round(mean(plot_WT_table),2));
txt_scn1a = 'Scn1a+/- Mean Corr: ' + string(round(mean(plot_scn1a_table),2));
text((mean(plot_WT_table)-0.05),(max(N)+15),txt_WT,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_scn1a_table)-0.05),(max(N)+5),txt_scn1a,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('Shuffled Neural Activity and Pupil','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
% Kolmogorov-Smirnov Test for Significant Difference Between Distributions

disp('Two Sample K-S Testing')
[h, p] = kstest2(plot_WT_table, plot_scn1a_table);
disp(h)
disp(p)
    
%% Plot Lag Histograms of all WT & Scn1a cells: Neural Activity and Pupil

window = 300;

% WT Neural Activity and Running
N = hist(WT_NP_lags, [-window:10:window]); % Create histogram of correlation coeffs
pd = fitdist(WT_NP_lags.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(WT_NP_lags)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-window:10:window],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.5;

hold on;

% Scn1a Neural Activity and Running
N = hist(scn1a_NP_lags,[-window:10:window]); % Create histogram of correlation coeffs
pd = fitdist(scn1a_NP_lags.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(scn1a_NR_lags)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-window:10:window],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.5;

hold on;
plot([mean(WT_NP_lags) mean(WT_NP_lags)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(scn1a_NP_lags) mean(scn1a_NP_lags)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_WT = 'WT Mean Lag: ' + string(round(mean(WT_NP_lags)/fs,2)) + ' s';
txt_scn1a = 'Scn1a+/- Mean Lag: ' + string(round(mean(scn1a_NP_lags)/fs,2)) + ' s';
text((mean(WT_NP_lags)-0.05),(max(N)),txt_WT,'HorizontalAlignment', 'right');
hold on;
text((mean(scn1a_NP_lags)-0.05),(max(N)-5),txt_scn1a,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Lag (s)','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('Lags: Neural Activity and Pupil','fontsize',14);
    xticks([-window:int16(fs)*2:window]);
    xticklabels([-window:int16(fs)*2:window]/fs);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
%% Plot Shuffled Lag Histograms of all WT & Scn1a cells: Neural Activity and Pupil

window = 300;

% WT Neural Activity and Running

plot_WT_table = WT_NP_lags_shuffle;

N = hist(plot_WT_table, [-window:10:window]); % Create histogram of correlation coeffs
pd = fitdist(plot_WT_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_WT_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-window:10:window],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.5;

hold on;

% Scn1a Neural Activity and Running

plot_scn1a_table = scn1a_NP_lags_shuffle;

N = hist(plot_scn1a_table,[-window:10:window]); % Create histogram of correlation coeffs
pd = fitdist(plot_scn1a_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_scn1a_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-window:10:window],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.5;

hold on;
plot([mean(plot_WT_table) mean(plot_WT_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(plot_scn1a_table) mean(plot_scn1a_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_WT = 'WT Mean Lag: ' + string(round(mean(plot_WT_table)/fs,2)) + ' s';
txt_scn1a = 'Scn1a+/- Mean Lag: ' + string(round(mean(plot_scn1a_table)/fs,2)) + ' s';
text((mean(plot_WT_table)-0.05),(max(N)+15),txt_WT,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_scn1a_table)-0.05),(max(N)+5),txt_scn1a,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Lag (s)','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('Shuffled Lags: Neural Activity and Pupil','fontsize',14);
    xticks([-window:int16(fs)*2:window]);
    xticklabels([-window:int16(fs)*2:window]/fs);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');    

%% WT Real and Shuffled Max Correlation: Neural Activity and Locomotion

% WT Real Data
plot_real_table = WT_NR;

N = hist(plot_real_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_real_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_real_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.6;

hold on;

% Scn1a Neural Activity and Running

plot_shuffle_table = WT_NR_shuffle;

N = hist(plot_shuffle_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_shuffle_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_shuffle_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b2.FaceAlpha = 0.3;

hold on;
plot([mean(plot_real_table) mean(plot_real_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(plot_shuffle_table) mean(plot_shuffle_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
txt_real = 'Real Mean Corr: ' + string(round(mean(plot_real_table),2));
txt_shuffle = 'Shuffled Mean Corr: ' + string(round(mean(plot_shuffle_table),2));
text((mean(plot_real_table)-0.05),(max(N)+15),txt_real,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_shuffle_table)-0.05),(max(N)+5),txt_shuffle,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('WT Neural Activity and Locomotion','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
% Kolmogorov-Smirnov Test for Significant Difference Between Distributions

disp('Two Sample K-S Testing')
[h, p] = kstest2(plot_real_table, plot_shuffle_table);
disp(h)
disp(p)
    
%% Scn1a Real and Shuffled Max Correlation: Neural Activity and Locomotion

% Scn1a Real Data
plot_real_table = scn1a_NR;

N = hist(plot_real_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_real_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_real_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b1.FaceAlpha = 0.6;

hold on;

% Scn1a Shuffled Data

plot_shuffle_table = scn1a_NR_shuffle;

N = hist(plot_shuffle_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_shuffle_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_shuffle_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.3;

hold on;
plot([mean(plot_real_table) mean(plot_real_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
plot([mean(plot_shuffle_table) mean(plot_shuffle_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_real = 'Real Mean Corr: ' + string(round(mean(plot_real_table),2));
txt_shuffle = 'Shuffled Mean Corr: ' + string(round(mean(plot_shuffle_table),2));
text((mean(plot_real_table)-0.05),(max(N)+15),txt_real,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_shuffle_table)-0.05),(max(N)+5),txt_shuffle,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('Scn1a Neural Activity and Locomotion','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
% Kolmogorov-Smirnov Test for Significant Difference Between Distributions

disp('Two Sample K-S Testing')
[h, p] = kstest2(plot_real_table, plot_shuffle_table);
disp(h)
disp(p)
    
    
%% WT Real and Shuffled Max Correlation: Neural Activity and Pupil

% WT Real Data
plot_real_table = WT_NP;

N = hist(plot_real_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_real_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_real_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.6;

hold on;

% Scn1a Shuffled Data

plot_shuffle_table = WT_NP_shuffle;

N = hist(plot_shuffle_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_shuffle_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_shuffle_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b2.FaceAlpha = 0.3;

hold on;
plot([mean(plot_real_table) mean(plot_real_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(plot_shuffle_table) mean(plot_shuffle_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
txt_real = 'Real Mean Corr: ' + string(round(mean(plot_real_table),2));
txt_shuffle = 'Shuffled Mean Corr: ' + string(round(mean(plot_shuffle_table),2));
text((mean(plot_real_table)-0.05),(max(N)+15),txt_real,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_shuffle_table)-0.05),(max(N)+5),txt_shuffle,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('WT Neural Activity and Pupil','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
% Kolmogorov-Smirnov Test for Significant Difference Between Distributions

disp('Two Sample K-S Testing')
[h, p] = kstest2(plot_real_table, plot_shuffle_table);
disp(h)
disp(p)
    
    
%% Scn1a Real and Shuffled Max Correlation: Neural Activity and Pupil

% Scn1a Real Data
plot_real_table = scn1a_NP;

N = hist(plot_real_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_real_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_real_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b1.FaceAlpha = 0.6;

hold on;

% Scn1a Shuffled Data

plot_shuffle_table = scn1a_NP_shuffle;

N = hist(plot_shuffle_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_shuffle_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_shuffle_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.3;

hold on;
plot([mean(plot_real_table) mean(plot_real_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
plot([mean(plot_shuffle_table) mean(plot_shuffle_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_real = 'Real Mean Corr: ' + string(round(mean(plot_real_table),2));
txt_shuffle = 'Shuffled Mean Corr: ' + string(round(mean(plot_shuffle_table),2));
text((mean(plot_real_table)-0.05),(max(N)+15),txt_real,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_shuffle_table)-0.05),(max(N)+5),txt_shuffle,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('Scn1a Neural Activity and Pupil','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
% Kolmogorov-Smirnov Test for Significant Difference Between Distributions

disp('Two Sample K-S Testing')
[h, p] = kstest2(plot_real_table, plot_shuffle_table);
disp(h)
disp(p)
    
%% Plot Zero Correlation Histograms of all WT & Scn1a cells: Neural Activity and Locomotion

% WT Neural Activity and Running
N = hist(WT_zero_NR, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(WT_zero_NR.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(WT_zero_NR)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.5;

hold on;

% Scn1a Neural Activity and Running
N = hist(scn1a_zero_NR,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(scn1a_zero_NR.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(scn1a_zero_NR)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.5;

hold on;
plot([mean(WT_zero_NR) mean(WT_zero_NR)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(scn1a_zero_NR) mean(scn1a_zero_NR)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_WT = 'WT mean Corr: ' + string(round(mean(WT_zero_NR),2));
txt_scn1a = 'Scn1a+/- mean Corr: ' + string(round(mean(scn1a_zero_NR),2));
text((mean(WT_zero_NR)-0.05),(max(N)),txt_WT,'HorizontalAlignment', 'right');
hold on;
text((mean(scn1a_zero_NR)-0.05),(max(N)-5),txt_scn1a,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('Neural Activity and Locomotion (Zero Lag)','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
% Kolmogorov-Smirnov Test for Significant Difference Between Distributions

disp('Two Sample K-S Testing')
[h, p] = kstest2(WT_zero_NR, scn1a_zero_NR);
disp(h)
disp(p)
    
disp('Wilcoxin Rank Sum Testing')
[h, p, stats] = ranksum(WT_zero_NR, scn1a_zero_NR);
disp(h)
disp(p)



%% WT Real and Shuffled Zero Correlation: Neural Activity and Locomotion

% WT Real Data
plot_real_table = WT_zero_NR;

N = hist(plot_real_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_real_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_real_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.6;

hold on;

% WT Shuffled Data

plot_shuffle_table = WT_zero_NR_shuffle;

N = hist(plot_shuffle_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_shuffle_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_shuffle_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b2.FaceAlpha = 0.3;

hold on;
plot([mean(plot_real_table) mean(plot_real_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(plot_shuffle_table) mean(plot_shuffle_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
txt_real = 'Real mean Corr: ' + string(round(mean(plot_real_table),2));
txt_shuffle = 'Shuffled mean Corr: ' + string(round(mean(plot_shuffle_table),2));
text((mean(plot_real_table)-0.05),(max(N)+15),txt_real,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_shuffle_table)-0.05),(max(N)+5),txt_shuffle,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('WT Neural Activity and Locomotion (Zero Lag)','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
% Test for Significant Difference Between Distributions

disp('Two Sample K-S Testing')
[h, p] = kstest2(plot_real_table, plot_shuffle_table);
disp(h)
disp(p)

disp('Wilcoxin Rank Sum Testing')
[h, p] = ranksum(plot_real_table, plot_shuffle_table);
disp(h)
disp(p)

% Find 5th and 95th percentile of shuffled correlations between pupil
% diameter and locomotion

pct5_nr_shuff = prctile(plot_shuffle_table, 5);
pct95_nr_shuff = prctile(plot_shuffle_table, 95);
    
%% Scn1a Real and Shuffled Zero Correlation: Neural Activity and Locomotion

% Scn1a Real Data
plot_real_table = scn1a_zero_NR;

N = hist(plot_real_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_real_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_real_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b1.FaceAlpha = 0.6;

hold on;

% Scn1a Shuffled Data

plot_shuffle_table = scn1a_zero_NR_shuffle;

N = hist(plot_shuffle_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_shuffle_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_shuffle_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.3;

hold on;
plot([mean(plot_real_table) mean(plot_real_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
plot([mean(plot_shuffle_table) mean(plot_shuffle_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_real = 'Real Mean Corr: ' + string(round(mean(plot_real_table),2));
txt_shuffle = 'Shuffled Mean Corr: ' + string(round(mean(plot_shuffle_table),2));
text((mean(plot_real_table)-0.05),(max(N)+15),txt_real,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_shuffle_table)-0.05),(max(N)+5),txt_shuffle,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('Scn1a Neural Activity and Locomotion (Zero Lag)','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
    
%% Plot Zero Correlation Histograms of all WT & Scn1a cells: Neural Activity and Pupil

% WT Neural Activity and Pupil
N = hist(WT_zero_NP, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(WT_zero_NP.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(WT_zero_NP)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.5;

hold on;

% Scn1a Neural Activity and Running
N = hist(scn1a_zero_NP,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(scn1a_zero_NP.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(scn1a_zero_NP)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.5;

hold on;
plot([mean(WT_zero_NP) mean(WT_zero_NP)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(scn1a_zero_NP) mean(scn1a_zero_NP)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_WT = 'WT mean Corr: ' + string(round(mean(WT_zero_NP),2));
txt_scn1a = 'Scn1a+/- mean Corr: ' + string(round(mean(scn1a_zero_NP),2));
text((mean(WT_zero_NP)-0.05),(max(N)),txt_WT,'HorizontalAlignment', 'right');
hold on;
text((mean(scn1a_zero_NP)-0.05),(max(N)-5),txt_scn1a,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('Neural Activity and Pupil (Zero Lag)','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');    
    
% Kolmogorov-Smirnov Test for Significant Difference Between Distributions

disp('Two Sample K-S Testing')
[h, p] = kstest2(WT_zero_NP, scn1a_zero_NP);
disp(h)
disp(p)

disp('Wilcoxin Rank Sum Testing')
[h, p] = ranksum(WT_zero_NP, scn1a_zero_NP);
disp(h)
disp(p)
    
%% WT Real and Shuffled Zero Correlation: Neural Activity and Pupil

% WT Real Data
plot_real_table = WT_zero_NP;

N = hist(plot_real_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_real_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_real_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.6;

hold on;

% WT Shuffled Data

plot_shuffle_table = WT_zero_NP_shuffle;

N = hist(plot_shuffle_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_shuffle_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_shuffle_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b2.FaceAlpha = 0.3;

hold on;
plot([mean(plot_real_table) mean(plot_real_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(plot_shuffle_table) mean(plot_shuffle_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
txt_real = 'Real Mean Corr: ' + string(round(mean(plot_real_table),2));
txt_shuffle = 'Shuffled Mean Corr: ' + string(round(mean(plot_shuffle_table),2));
text((mean(plot_real_table)-0.05),(max(N)+15),txt_real,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_shuffle_table)-0.05),(max(N)+5),txt_shuffle,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('WT Neural Activity and Pupil (Zero Lag)','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
% Test for Significant Difference Between Distributions

disp('Two Sample K-S Testing')
[h, p] = kstest2(plot_real_table, plot_shuffle_table);
disp(h)
disp(p)

disp('Wilcoxin Rank Sum Testing')
[h, p] = ranksum(plot_real_table, plot_shuffle_table);
disp(h)
disp(p)
    
%% Scn1a Real and Shuffled Zero Correlation: Neural Activity and Pupil

% Scn1a Real Data
plot_real_table = scn1a_zero_NP;

N = hist(plot_real_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_real_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_real_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b1.FaceAlpha = 0.6;

hold on;

% Scn1a Shuffled Data

plot_shuffle_table = scn1a_zero_NP_shuffle;

N = hist(plot_shuffle_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_shuffle_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_shuffle_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.3;

hold on;
plot([mean(plot_real_table) mean(plot_real_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
plot([mean(plot_shuffle_table) mean(plot_shuffle_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_real = 'Real Mean Corr: ' + string(round(mean(plot_real_table),2));
txt_shuffle = 'Shuffled Mean Corr: ' + string(round(mean(plot_shuffle_table),2));
text((mean(plot_real_table)-0.05),(max(N)+15),txt_real,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_shuffle_table)-0.05),(max(N)+5),txt_shuffle,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('Scn1a Neural Activity and Pupil (Zero Lag)','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');    
    
%% Mask Out Running Bouts from dF and pupil data

speed_th = 10; % Set speed threshold to count as locomotion bout
win_l = 2; % Set window for moving mean smoothing by Inst_Speed_Mask fxn (sec)
win_count = int16(10*fs); % Set minimum length of quiet bout

for i = 1:length(data)
    
    % Use Inst_Spped_Mask function to get an array with 1 when mouse is
    % quiet and 0 when mouse is running
    [mask_loc]=Inst_Speed_Mask(data(i).speed,data(i).fs,speed_th,win_l);
    %figure;plot(inst_speed_r);
    %hold on;plot(mask_loc*speed_th);
    
    [wid_p, initcr_p, fincr_p, midcr_p] = pulsewidth(double(not(mask_loc)));
    % wid_p = width of each running bout
    % initcr_p = start index of each running bout
    % fincr_p = end index of each running bout
    
    mask_mod = zeros(size(mask_loc));
    ct = 1; % counter for running bouts
    st_pt = []; % Array of starting indices for running bouts
    dF_quiet = {}; % Array of dF/F values during quiesence (non-running)
    pupil_quiet = {}; % Array of pupil values during quiesence
    
    
    
    for jj=1:length(wid_p) % looping through each running bout
        if jj == 1
            % make sure that you also get the beginning of the recoring if
            % the mouse is quiet at the beginning
            if initcr_p(jj) > 2
                
                % first quiet bout
                st_pt(ct) = 1; % first chunk
                start_index = 1;
                fin_index = int16(initcr_p(1)); % set end index as start of 1st running bout
                if (fin_index - start_index) > win_count
                    dF_quiet{ct} = data(i).dF(:,1:fin_index); % assign dF values during 1st quiet bout
                    pupil_quiet{ct} = data(i).pupil(1:fin_index); % assign pupil values during 1st quiet
                    ct = ct + 1;
                end
                
                % second quiet bout
                st_pt(ct) = int16(fincr_p(1)); % set start index as end of 1st running bout
                start_index = int16(fincr_p(jj));
                fin_index = int16(initcr_p(jj+1)); % set end index as start of 2nd running bout
                if (fin_index - start_index) > win_count
                    dF_quiet{ct} = data(i).dF(:,start_index:fin_index); % assign dF values during 1st quiet bout
                    pupil_quiet{ct} = data(i).pupil(start_index:fin_index); % assign pupil values during 1st quiet
                    ct = ct + 1;
                end
                
            else % mouse running at beginning of recording
                
                % first quiet bout
                st_pt(ct) = int16(fincr_p(1)); % set start index as end of 1st running bout
                start_index = int16(fincr_p(jj));
                fin_index = int16(initcr_p(jj+1)); % set end index as start of 2nd running bout
                if (fin_index - start_index) > win_count
                    dF_quiet{ct} = data(i).dF(:,start_index:fin_index); % assign dF values during 1st quiet bout
                    pupil_quiet{ct} = data(i).pupil(start_index:fin_index); % assign pupil values during 1st quiet
                    ct = ct + 1;
                end
                
            end
            
        elseif initcr_p(jj) == max(initcr_p)  % if last running bout in recording
            
            % last quiet bout
            st_pt(ct) = int16(fincr_p(jj)); % set start index as end of 1st running bout
            start_index = int16(fincr_p(jj));
            fin_index = width(data(i).dF); % set end index as start of 2nd running bout
            if (fin_index - start_index) > win_count
                dF_quiet{ct} = data(i).dF(:,start_index:fin_index); % assign dF values during 1st quiet bout
                pupil_quiet{ct} = data(i).pupil(start_index:fin_index); % assign pupil values during 1st quiet
            end
            
        else
            
            % middle quiet bout
            st_pt(ct) = int16(fincr_p(jj)); % set start index as end of 1st running bout
            start_index = int16(fincr_p(jj));
            fin_index = int16(initcr_p(jj+1)); % set end index as start of 2nd running bout
            if (fin_index - start_index) > win_count
                dF_quiet{ct} = data(i).dF(:,start_index:fin_index); % assign dF values during 1st quiet bout
                pupil_quiet{ct} = data(i).pupil(start_index:fin_index); % assign pupil values during 1st quiet
                ct = ct + 1;
            end
            
        end
    end
    
    % dF array Structure: (locomotion bout, cell number, time)
    data(i).dF_quiet = dF_quiet; % Add to data table
    data(i).pupil_quiet = pupil_quiet;
    
end

%% Calculate cross correlation between neural activity and pupil during quiescence only

for i = 1:length(data) % loop over each recording
    
    clear corr_NP_avg lags_NP_avg corr_quiet_avg corr_quiet_shuffle_avg
            
    for ii = 1:size(data(i).dF_quiet,2); % loop over each quiet bout
        
        clear corr_np_quiet lags_np_quiet...
        corr_NP lags_NP corrs_quiet
    
        pupil = zscore(fillmissing(data(i).pupil_quiet{ii},'linear'));
        neural_array = data(i).dF_quiet{ii}; % all cells during quiet bout
    
        for iii=1:size(neural_array,1) % loop over each cell
        
            % Data files for each cell for alignment
            neural = zscore(neural_array(iii,:));
            
            % calculate cross correlation of neural activity and pupil
            [corr_np_quiet, lags_np_quiet] = xcorr(neural, pupil, 300, 'normalized');
            
            % Get max correlation and lag for that max correlation value and
            % assign to table
            [c_NP, l_NP] = max(abs(corr_np_quiet));
            corr_NP(iii) = corr_np_quiet(l_NP); % get max correlation value
            lags_NP(iii) = lags_np_quiet(l_NP); % get lag (index) of max correlation
       
        end
        
        corr_NP_avg(ii,:,:) = corr_NP;
        lags_NP_avg(ii,:,:) = lags_NP;
        
        % Calculate cell-cell correlations during quiet
        [rho, pval] = corr(zscore(data(i).dF_quiet{ii},1,2)', zscore(data(i).dF_quiet{ii},1,2)');
        corr_quiet_avg(ii,:,:) = rho;
        [rho, pval] = corr(zscore(data(i).dF_quiet{ii},1,2)', zscore(flip(data(i).dF_quiet{ii},2),1,2)');
        corr_quiet_shuffle_avg(ii,:,:) = rho;
    
    end
    
    % Add average correlation data for each cell across quiet bouts
    data(i).corr_NP_quiet = squeeze(mean(corr_NP_avg,1));
    data(i).lags_NP_quiet = squeeze(mean(lags_NP_avg,1));
    data(i).corrs_quiet = squeeze(mean(corr_quiet_avg,1));
    data(i).corrs_quiet_shuffle = squeeze(mean(corr_quiet_shuffle_avg,1));
    
end

%% Pupil Correlation Histograms by Genotype during Quiet

% Create tables with correlation coefficients for all WT/Scn1a mice

WT_NP_quiet = [];
WT_NP_quiet_lags = [];
scn1a_NP_quiet = [];
scn1a_NP_quiet_lags = [];


for i=1:length(data) % loop over each recording
        
    if ismember(data(i).mouse(1),WT)
        WT_NP_quiet = [WT_NP_quiet; data(i).corr_NP_quiet];
        WT_NP_lags = [WT_NP_quiet_lags; data(i).lags_NP_quiet];
    elseif ismember(data(i).mouse(1),scn1a)
        scn1a_NP_quiet = [scn1a_NP_quiet; data(i).corr_NP_quiet];
        scn1a_NP_quiet_lags = [scn1a_NP_quiet_lags; data(i).lags_NP_quiet];

    end
    
end

%% Plot Correlation Histograms of all WT & Scn1a cells: Neural Activity and Pupil during Quiet Only

% WT Neural Activity and Running
N = hist(WT_NP_quiet, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(WT_NP_quiet,'normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(WT_NP_quiet)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.5;

hold on;

% Scn1a Neural Activity and Running
N = hist(scn1a_NP_quiet,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(scn1a_NP_quiet,'normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(scn1a_NP_quiet)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.5;

hold on;
plot([mean(WT_NP_quiet) mean(WT_NP_quiet)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(scn1a_NP_quiet) mean(scn1a_NP_quiet)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_WT = 'WT Mean Corr: ' + string(round(mean(WT_NP_quiet),2));
txt_scn1a = 'Scn1a+/- Mean Corr: ' + string(round(mean(scn1a_NP_quiet),2));
text((mean(WT_NP_quiet)-0.05),(max(N)),txt_WT,'HorizontalAlignment', 'right');
hold on;
text((mean(scn1a_NP_quiet)-0.05),(max(N)-5),txt_scn1a,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cells','fontsize',14);
    title('Neural Activity and Pupil (Quiet)','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');

%% Calculate dF/F Change at Locomotion Onset

speed_th = 10; % Set speed threshold to count as locomotion bout
win_l = 2; % Set window for moving mean smoothing by Inst_Speed_Mask fxn (sec)
win1 = 4; % Set window for time before running onset (sec)
win2 = 6; % Set window for time after running onset (sec)

for i = 1:length(data)
    
    % Use Inst_Spped_Mask function to get an array with 1 when mouse is
    % quiet and 0 when mouse is running
    [mask_loc]=Inst_Speed_Mask(data(i).speed,data(i).fs,speed_th,win_l);
    %figure;plot(inst_speed_r);
    %hold on;plot(mask_loc*speed_th);
    
    [wid_p, initcr_p, fincr_p, midcr_p] = pulsewidth(double(not(mask_loc)));
    % wid_p = width of each running bout
    % initcr_p = start index of each running bout
    % fincr_p = end index of each running bout
    
    mask_mod = zeros(size(mask_loc));
    dF_trans = []; % Array of dF/F values during onset of locomotion windows (all bouts)
    speed_trans = []; % Array of speed values during onset of locomotion windows
    ct = 0; % counter for running bouts
    st_pt = []; % Array of starting indices for running bouts
    for jj=1:length(wid_p) % looping through each running bout
        % if start of recording or greater than 1 time window has passed,
        % count as new running bout
        if (jj==1 || (round(initcr_p(jj))-round(fincr_p(jj-1)))>(win1+win2)*fs) %round(initcr_p(jj-1)(round(initcr_p(jj))-st_pt(ct))>(win1+win2)*fs
            ct=ct+1;
        st_pt(ct) = round(initcr_p(jj));
        speed_trans(ct,:) = data(i).speed(st_pt(ct)-win1*fs+1:st_pt(ct)+win2*fs)-mean(data(i).speed(st_pt(ct)-win1*fs+1:st_pt(ct)));
        dF_trans(ct,:,:) = (data(i).dF(:,st_pt(ct)-win1*fs+1:st_pt(ct)+win2*fs)...
            -repmat(mean(data(i).dF(:,st_pt(ct)-win1*fs+1:st_pt(ct)),2),1,size(st_pt(ct)-win1*fs+1:st_pt(ct)+win2*fs,2)))*100;
        mask_mod(st_pt(ct)-win1*fs+1:st_pt(ct)+win2*fs)=1;
        end
    end
    
    % dF array Structure: (locomotion bout, cell number, time)
    data(i).speed_trans = speed_trans; % Add to data table
    data(i).dF_trans = dF_trans;
    
end

% Create arrays with all WT and Scn1a+/- cell data at onset of locomotion

WT_trans = []; % Store WT dF/F data at onset of locomotion (by each locomotion
% bout for each cell)
WT_transm = [];  % Store WT dF/F data at onset of locomotion (mean of locomotion
% bouts for each cell)
scn1a_trans = []; % Store Scn1a dF/F data at onset of locomotion 
% bout for each cell)
scn1a_transm = [];  % Store Scn1a dF/F data at onset of locomotion (mean of locomotion
% bouts for each cell)

for i=1:length(data)
    
    if ismember(data(i).mouse(1),WT)
        for j = 1:size(data(i).dF_trans,1) % Loop over each locomotion bout
            WT_trans = [WT_trans; squeeze(data(i).dF_trans(j,:,:))]; % Add each cell to WT_trans array
        end
        WT_transm = [WT_transm; squeeze(mean(data(i).dF_trans,1))];
    elseif ismember(data(i).mouse(1),scn1a)
        for j = 1:size(data(i).dF_trans,1) % Loop over each locomotion bout
            scn1a_trans = [scn1a_trans; squeeze(data(i).dF_trans(j,:,:))]; % Add each cell to WT_trans array
        end
        scn1a_transm = [scn1a_transm; squeeze(mean(data(i).dF_trans,1))];
    end
    
end

%% Calculate 95% Confidence Intervals for dF/F at Onset of Locomotion
% Counting all locomotion bouts as individual n

nBoot = 2000; % Number of bootstraps for resampling

[WT_bci,WT_bmeans] = bootci(nBoot,{@mean,WT_trans},'alpha',0.05,'type','per');
[scn1a_bci,scn1a_bmeans] = bootci(nBoot,{@mean,scn1a_trans},'alpha',0.05,'type','per');

%% Calculate 95% Confidence Intervals for dF/F at Onset of Locomotion
% Counting each cell as individual n

nBoot = 2000; % Number of bootstraps for resampling

[WT_bcim,WT_bmeansm] = bootci(nBoot,{@mean,WT_transm},'alpha',0.05,'type','per');
[scn1a_bcim,scn1a_bmeansm] = bootci(nBoot,{@mean,scn1a_transm},'alpha',0.05,'type','per');


%% Plot Mean dF/F at Onset of Locomotion
% Counting all locomotion bouts as individual n

set(gcf,'position',[100,100,540,400]);
set(gca, 'FontName', 'Arial');

time_trans = (-win1*fs+1:win2*fs)./fs;

% Plot Scn1a mean dF/F
plot(time_trans,mean(scn1a_trans,1),'Color',scn1a_color,'LineWidth',4); hold on;
% Plot SEM as shading around line
% curve1 = mean(scn1a_trans,1) + std(scn1a_trans,1)./sqrt(height(scn1a_trans));
% curve2 = mean(scn1a_trans,1) - std(scn1a_trans,1)./sqrt(height(scn1a_trans));
% Plot 95% CI as shading around line
curve1 = scn1a_bci(1,:);
curve2 = scn1a_bci(2,:);
x2 = [time_trans, fliplr(time_trans)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween,scn1a_color,'FaceAlpha',0.3,'EdgeAlpha',0);
hold on;

% Plot WT mean dF/F
plot(time_trans,mean(WT_trans,1),'k','LineWidth',4); hold on;
% Plot SEM as shading around line
curve1 = WT_bci(1,:);
curve2 = WT_bci(2,:);
x2 = [time_trans, fliplr(time_trans)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'k','FaceAlpha',0.3,'EdgeAlpha',0);
hold on;

title('Mean dF/F at Onset of Locomotion');
xlabel('Time (sec)');
ylabel('dF/F0 (%)');

ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');

%% Plot Mean dF/F at Onset of Locomotion
% Counting each cell as individual n

set(gcf,'position',[100,100,540,400]);
set(gca, 'FontName', 'Arial');

time_trans = (-win1*fs+1:win2*fs)./fs;

% Plot Scn1a mean dF/F
plot(time_trans,mean(scn1a_transm,1),'Color',scn1a_color,'LineWidth',4); hold on;
% Plot SEM as shading around line
% curve1 = mean(scn1a_trans,1) + std(scn1a_trans,1)./sqrt(height(scn1a_trans));
% curve2 = mean(scn1a_trans,1) - std(scn1a_trans,1)./sqrt(height(scn1a_trans));
% Plot 95% CI as shading around line
curve1 = scn1a_bcim(1,:);
curve2 = scn1a_bcim(2,:);
x2 = [time_trans, fliplr(time_trans)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween,scn1a_color,'FaceAlpha',0.3,'EdgeAlpha',0);
hold on;

% Plot WT mean dF/F
plot(time_trans,mean(WT_transm,1),'k','LineWidth',4); hold on;
% Plot SEM as shading around line
curve1 = WT_bcim(1,:);
curve2 = WT_bcim(2,:);
x2 = [time_trans, fliplr(time_trans)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'k','FaceAlpha',0.3,'EdgeAlpha',0);
hold on;

title('Mean dF/F at Onset of Locomotion');
xlabel('Time (sec)');
ylabel('dF/F0 (%)');
%% Create sorted arrays of all cells at locomotion onset

WT_norm = WT_transm;
scn1a_norm = scn1a_transm;

% Create arrays to sort by diff in dF/F 1 sec before and after onset of
% locomotion
% WT_diffs = mean(WT_norm(:,(win1-1)*fs:win1*fs),2) - mean(WT_norm(:,win1*fs:(win1+0.25)*fs),2);
% scn1a_diffs = mean(scn1a_norm(:,(win1-1)*fs:win1*fs),2) - mean(scn1a_norm(:,win1*fs:(win1+0.25)*fs),2);
WT_diffs = mean(WT_norm(:,win1*fs+1:end),2);
scn1a_diffs = mean(scn1a_norm(:,win1*fs+1:end),2);

% Sort locomotion diffs in descending order and store sorting index
%[WT_diffs,WT_idx] = sort(WT_diffs,'descend');
%[scn1a_diffs,scn1a_idx] = sort(scn1a_diffs,'descend');

WT_norm = horzcat(WT_diffs,WT_norm);
scn1a_norm = horzcat(scn1a_diffs,scn1a_norm);
WT_norm = sortrows(WT_norm, 1, 'descend');
scn1a_norm = sortrows(scn1a_norm, 1, 'descend');
WT_norm = WT_norm(:,2:end);
scn1a_norm = scn1a_norm(:,2:end);

%% Plot colormesh of all WT cells at locomotion onset

set(gcf,'position',[100,100,540,400]);
set(gca, 'FontName', 'Arial');

h = heatmap(WT_norm, 'Colormap', turbo,...
    'Colorlimits',[prctile(WT_norm(:),15), prctile(WT_norm(:),85)]);
h.GridVisible = 'off';

xidx = 1:width(WT_norm);
xidx = rem(xidx,25)==0;

yidx = 1:height(WT_norm);
yidx = rem(yidx,100)==0;

h.XDisplayLabels(~xidx) = {''};
h.YDisplayLabels(~yidx) = {''};

xlabel('Time')
ylabel('Cell Number')
title('WT dF/F at Locomotion Onset')

%% Plot colormesh of all Scn1a cells at locomotion onset

set(gcf,'position',[100,100,540,400]);
set(gca, 'FontName', 'Arial');

figure;
h = heatmap(scn1a_norm, 'Colormap', turbo,...
    'Colorlimits',[prctile(WT_norm(:),15), prctile(WT_norm(:),85)]);
h.GridVisible = 'off';

xidx = 1:width(scn1a_norm);
xidx = rem(xidx,25)==0;

yidx = 1:height(scn1a_norm);
yidx = rem(yidx,100)==0;

h.XDisplayLabels(~xidx) = {''};
h.YDisplayLabels(~yidx) = {''};

xlabel('Time')
ylabel('Cell Number')
title('Scn1a+/- dF/F at Locomotion Onset')

%% Calculate dF/F Change at Locomotion Onset for Responders

% Create arrays with all WT and Scn1a+/- cell data at onset of locomotion
% Positive responders
WT_trans_resp = []; % Store WT dF/F data at onset of locomotion (by each locomotion
% bout for each cell)
WT_transm_resp = [];  % Store WT dF/F data at onset of locomotion (mean of locomotion
% bouts for each cell)
scn1a_trans_resp = []; % Store Scn1a dF/F data at onset of locomotion 
% bout for each cell)
scn1a_transm_resp = [];  % Store Scn1a dF/F data at onset of locomotion (mean of locomotion
% bouts for each cell)

% Counters for numbers of each type of responder
total_WT_responders = 0;
total_WT_neg_responders = 0;
total_WT_non_responders = 0;
total_scn1a_responders = 0;
total_scn1a_neg_responders = 0;
total_scn1a_non_responders = 0;

for i=1:length(data)
    responders = mean(data(i).zero_corr_NR,1) > pct95_nr_shuff;
    neg_responders = mean(data(i).zero_corr_NR,1) < pct5_nr_shuff;
    non_responders = mean(data(i).zero_corr_NR,1) > pct5_nr_shuff & mean(data(i).zero_corr_NR,1) < pct95_nr_shuff;
    
    if ismember(data(i).mouse(1),WT)
        for j = 1:size(data(i).dF_trans,1) % Loop over each locomotion bout           
            WT_trans_resp = [WT_trans_resp; squeeze(data(i).dF_trans(j,responders,:))]; % Add each cell to WT_trans array           
        end
        WT_transm_resp = [WT_transm_resp; squeeze(mean(data(i).dF_trans(:,responders,:),1))];
        
        % Count total number of responders, negative, and non-responders
        total_WT_responders = total_WT_responders + sum(responders == 1);
        total_WT_neg_responders = total_WT_neg_responders + sum(neg_responders == 1);
        total_WT_non_responders = total_WT_non_responders + sum(non_responders == 1);
        
    elseif ismember(data(i).mouse(1),scn1a)
        for j = 1:size(data(i).dF_trans,1) % Loop over each locomotion bout
            scn1a_trans_resp = [scn1a_trans_resp; squeeze(data(i).dF_trans(j,responders,:))]; % Add each cell to WT_trans array
        end
        scn1a_transm_resp = [scn1a_transm_resp; squeeze(mean(data(i).dF_trans(:,responders,:),1))];
        
        % Count total number of responders, negative, and non-responders
        total_scn1a_responders = total_scn1a_responders + sum(responders == 1);
        total_scn1a_neg_responders = total_scn1a_neg_responders + sum(neg_responders == 1);
        total_scn1a_non_responders = total_scn1a_non_responders + sum(non_responders == 1);
    end
    
    
end

%% Percent of Responders, Non-responders, Negative Responders

WT_total_cells = total_WT_responders + total_WT_non_responders + total_WT_neg_responders
scn1a_total_cells = total_scn1a_responders + total_scn1a_non_responders + total_scn1a_neg_responders

disp("WT Responders: " + string(total_WT_responders / WT_total_cells * 100) + "%")
disp("WT Non-Responders: " + string(total_WT_non_responders / WT_total_cells * 100) + "%")
disp("WT Negative Responders: " + string(total_WT_neg_responders / WT_total_cells * 100) + "%")

disp("Scn1a+/- Responders: " + string(total_scn1a_responders / scn1a_total_cells * 100) + "%")
disp("Scn1a+/- Non-Responders: " + string(total_scn1a_non_responders / scn1a_total_cells * 100) + "%")
disp("Scn1a+/- Negative Responders: " + string(total_scn1a_neg_responders / scn1a_total_cells * 100) + "%")

% Fisher's Exact Test for Responders
[h, p] = fishertest([total_WT_responders,(WT_total_cells-total_WT_responders);total_scn1a_responders,(scn1a_total_cells-total_scn1a_responders)]);

disp("Fisher's Exact Test for Responders");
disp(p)

[h, p] = fishertest([total_WT_neg_responders,(WT_total_cells-total_WT_neg_responders);total_scn1a_neg_responders,(scn1a_total_cells-total_scn1a_neg_responders)]);

disp("Fisher's Exact Test for Negative Responders");
disp(p)

%% Plot Mean dF/F at Onset of Locomotion for Responders
% Counting each cell as individual n

%%% Calculate 95% Confidence Intervals for dF/F at Onset of Locomotion for Responders
% Counting each cell as individual n

nBoot = 2000; % Number of bootstraps for resampling

[WT_bcim_resp,WT_bmeansm_resp] = bootci(nBoot,{@mean,WT_transm_resp},'alpha',0.05,'type','per');
[scn1a_bcim_resp,scn1a_bmeansm_resp] = bootci(nBoot,{@mean,scn1a_transm_resp},'alpha',0.05,'type','per');

set(gcf,'position',[100,100,540,400]);
set(gca, 'FontName', 'Arial');

time_trans = (-win1*fs+1:win2*fs)./fs;

% Plot Scn1a mean dF/F
plot(time_trans,mean(scn1a_transm_resp,1),'Color',scn1a_color,'LineWidth',4); hold on;
% Plot SEM as shading around line
% curve1 = mean(scn1a_trans,1) + std(scn1a_trans,1)./sqrt(height(scn1a_trans));
% curve2 = mean(scn1a_trans,1) - std(scn1a_trans,1)./sqrt(height(scn1a_trans));
% Plot 95% CI as shading around line
curve1 = scn1a_bcim_resp(1,:);
curve2 = scn1a_bcim_resp(2,:);
x2 = [time_trans, fliplr(time_trans)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween,scn1a_color,'FaceAlpha',0.3,'EdgeAlpha',0);
hold on;

% Plot WT mean dF/F
plot(time_trans,mean(WT_transm_resp,1),'k','LineWidth',4); hold on;
% Plot SEM as shading around line
curve1 = WT_bcim_resp(1,:);
curve2 = WT_bcim_resp(2,:);
x2 = [time_trans, fliplr(time_trans)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'k','FaceAlpha',0.3,'EdgeAlpha',0);
hold on;

title('Mean dF/F at Onset of Locomotion (Responders Only)');
xlabel('Time (sec)');
ylabel('dF/F0 (%)');

%% Correlation Between All Cells

% All cell in data frame for correlation matrix of all cells in recording
for i = 1:length(data)

    [rho, pval] = corr(zscore(data(i).dF,1,2)', zscore(data(i).dF,1,2)');
    data(i).corrs = rho;
    [rho, pval] = corr(zscore(data(i).dF,1,2)', zscore(flip(data(i).dF,2),1,2)');
    data(i).corrs_shuffle = rho;
    %heatmap(rho);

end

%% Clustering Based on Correlation Using Hierarchical Custering

% corrMat = rho - eye(size(rho)); % remove diagonal elements
% dissimilarity = 1 - corrMat(find(corrMat)); % convert to a vector (as pdist)
% 
% % Cutoff for Correlation
% cutoff = 0.5;

%save_dir_n = "C:\Users\liebergals\OneDrive - Children's Hospital of Philadelphia\Goldberg Lab\Sophie\figures\ndnf_scn1a\2P\arousal\ndnf_clustering_figs";

for i = 1:length(data)

    % Perform complete linkage clustering
    Z = linkage(data(i).corrs,'complete');

    % Group data into clusters
    groups = cluster(Z, 'maxclust', 2, 'criterion', 'distance');

    % Plot Sorted Based on Groups
    subplot(3,1,1)
    if size(data(i).dF(groups==1,:),1) > 1
        h = pcolor(data(i).dF(groups==1,:));
        set(h, 'EdgeColor', 'none');
        colormap(turbo);
        title('Cluster 1');
    end
    
    if size(data(i).dF(groups==2,:),1) > 1
        subplot(3,1,2)
        h = pcolor(data(i).dF(groups==2,:));
        set(h, 'EdgeColor', 'none');
        colormap(turbo);
        title('Cluster 2');
    end
    
%     if size(data(i).dF(groups==3,:),1) > 1
%         subplot(4,1,3)
%         h = pcolor(data(i).dF(groups==3,:));
%         set(h, 'EdgeColor', 'none');
%         colormap(turbo);
%         title('Cluster 3');
%     end
    
    subplot(3,1,3)
    plot(data(i).speed);
    
    filename = append(data(i).name, '_clusters.tif');
    
    %saveas(gcf,fullfile(save_dir_n,filename));

end


%% Plot Histogram of all Cell-Cell Correlations by Genotype

WT_pcorrs = [];
scn1a_pcorrs = [];

for i = 1:length(data)
    
    lower_corr = [];
    lower_corr = tril(data(i).corrs, -1);
    lower_corr = reshape(lower_corr,[],1);
    lower_corr(lower_corr == 0) = [];
    
    if ismember(data(i).mouse(1),WT)
        WT_pcorrs = [WT_pcorrs; lower_corr];
    elseif ismember(data(i).mouse(1),scn1a)
        scn1a_pcorrs = [scn1a_pcorrs; lower_corr];
    end
    
end

% % Create x vectors for swarm plots
% WT_x = ones(length(WT_pcorrs),1);
% scn1a_x = 2*ones(length(scn1a_pcorrs),1);
% 
% % Plot all correlation coefficients as swarm charts with overlaid box plots
% swarmchart(WT_x,WT_pcorrs,'k');
% hold on
% swarmchart(scn1a_x,scn1a_pcorrs,'r');
% hold on
% g1 = repmat({'WT'},length(WT_x),1);
% g2 = repmat({'Scn1a'},length(scn1a_x),1);
% g = [g1; g2];
% boxplot([WT_pcorrs;scn1a_pcorrs],g,'Colors','w');
% title('Correlation Coefficients Between Ndnf+ Cells');
% xlabel('Genotype');
% ylabel('Correlation Coefficient');
% 
% WT
plot_WT_table = WT_pcorrs';

N = hist(plot_WT_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_WT_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_WT_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.5;

hold on;

% Scn1a

plot_scn1a_table = scn1a_pcorrs';

N = hist(plot_scn1a_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_scn1a_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_scn1a_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.5;

hold on;
plot([mean(plot_WT_table) mean(plot_WT_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(plot_scn1a_table) mean(plot_scn1a_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_WT = 'WT Mean Corr: ' + string(round(mean(plot_WT_table),2));
txt_scn1a = 'Scn1a+/- Mean Corr: ' + string(round(mean(plot_scn1a_table),2));
text((mean(plot_WT_table)-0.05),(max(N)+25),txt_WT,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_scn1a_table)-0.05),(max(N)-50),txt_scn1a,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cell Pairs','fontsize',14);
    title('Correlation Between Pairs of Ndnf-INs','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
       
% Test for Significant Difference Between Distributions

disp('Two Sample K-S Testing')
[h, p] = kstest2(plot_WT_table, plot_scn1a_table);
disp(h)
disp(p)

disp('Wilcoxin Rank Sum Testing')
[h, p] = ranksum(plot_WT_table, plot_scn1a_table);
disp(h)
disp(p)

%% Plot Shuffled Histogram of all Cell-Cell Correlations by Genotype

WT_pcorrs_shuffle = [];
scn1a_pcorrs_shuffle = [];

for i = 1:length(data)
    
    lower_corr = [];
    lower_corr = tril(data(i).corrs_shuffle, -1);
    lower_corr = reshape(lower_corr,[],1);
    lower_corr(lower_corr == 0) = [];
    
    if ismember(data(i).mouse(1),WT)
        WT_pcorrs_shuffle = [WT_pcorrs_shuffle; lower_corr];
    elseif ismember(data(i).mouse(1),scn1a)
        scn1a_pcorrs_shuffle = [scn1a_pcorrs_shuffle; lower_corr];
    end
    
end

% % Create x vectors for swarm plots
% WT_x = ones(length(WT_pcorrs),1);
% scn1a_x = 2*ones(length(scn1a_pcorrs),1);
% 
% % Plot all correlation coefficients as swarm charts with overlaid box plots
% swarmchart(WT_x,WT_pcorrs,'k');
% hold on
% swarmchart(scn1a_x,scn1a_pcorrs,'r');
% hold on
% g1 = repmat({'WT'},length(WT_x),1);
% g2 = repmat({'Scn1a'},length(scn1a_x),1);
% g = [g1; g2];
% boxplot([WT_pcorrs;scn1a_pcorrs],g,'Colors','w');
% title('Correlation Coefficients Between Ndnf+ Cells');
% xlabel('Genotype');
% ylabel('Correlation Coefficient');
% 
% WT Neural Activity and Running
plot_WT_table = WT_pcorrs_shuffle';

N = hist(plot_WT_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_WT_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_WT_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.5;

hold on;

% Scn1a Neural Activity and Running

plot_scn1a_table = scn1a_pcorrs_shuffle';

N = hist(plot_scn1a_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_scn1a_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_scn1a_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.5;

hold on;
plot([mean(plot_WT_table) mean(plot_WT_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(plot_scn1a_table) mean(plot_scn1a_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_WT = 'WT Mean Corr: ' + string(round(mean(plot_WT_table),2));
txt_scn1a = 'Scn1a+/- Mean Corr: ' + string(round(mean(plot_scn1a_table),2));
text((mean(plot_WT_table)-0.05),(max(N)+25),txt_WT,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_scn1a_table)-0.05),(max(N)-50),txt_scn1a,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cell Pairs','fontsize',14);
    title('Shuffled Correlation Between Pairs of Ndnf-INs','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
% Test for Significant Difference Between Distributions

disp('Two Sample K-S Testing')
[h, p] = kstest2(plot_WT_table, plot_scn1a_table);
disp(h)
disp(p)

disp('Wilcoxin Rank Sum Testing')
[h, p] = ranksum(plot_WT_table, plot_scn1a_table);
disp(h)
disp(p)
    
%% Plot WT Real and Shuffled Cell-Cell Correlations

% WT Real Data
plot_real_table = WT_pcorrs';

N = hist(plot_real_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_real_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_real_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.6;

hold on;

% Scn1a Neural Activity and Running

plot_shuffle_table = WT_pcorrs_shuffle';

N = hist(plot_shuffle_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_shuffle_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_shuffle_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b2.FaceAlpha = 0.3;

hold on;
plot([mean(plot_real_table) mean(plot_real_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(plot_shuffle_table) mean(plot_shuffle_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
txt_real = 'Real Mean Corr: ' + string(round(mean(plot_real_table),2));
txt_shuffle = 'Shuffled Mean Corr: ' + string(round(mean(plot_shuffle_table),2));
text((mean(plot_real_table)-0.05),(max(N)+15),txt_real,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_shuffle_table)-0.05),(max(N)-120),txt_shuffle,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cell-Cell Pairs','fontsize',14);
    title('WT Cell-Cell Correlations','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
    % Test for Significant Difference Between Distributions

disp('Two Sample K-S Testing')
[h, p] = kstest2(plot_real_table, plot_shuffle_table);
disp(h)
disp(p)

disp('Wilcoxin Rank Sum Testing')
 [h, p] = ranksum(plot_real_table, plot_shuffle_table);
 disp(h)
 disp(p)
    
%% Plot Scn1a Real and Shuffled Cell-Cell Correlations

% Scn1a Real Data
plot_real_table = scn1a_pcorrs';

N = hist(plot_real_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_real_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_real_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b1.FaceAlpha = 0.6;

hold on;

% Scn1a Shuffled Data

plot_shuffle_table = scn1a_pcorrs_shuffle';

N = hist(plot_shuffle_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_shuffle_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_shuffle_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.3;

hold on;
plot([mean(plot_real_table) mean(plot_real_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
plot([mean(plot_shuffle_table) mean(plot_shuffle_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_real = 'Real Mean Corr: ' + string(round(mean(plot_real_table),2));
txt_shuffle = 'Shuffled Mean Corr: ' + string(round(mean(plot_shuffle_table),2));
text((mean(plot_real_table)-0.05),(max(N)+15),txt_real,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_shuffle_table)-0.05),(max(N)-140),txt_shuffle,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cell-Cell Pairs','fontsize',14);
    title('Scn1a Cell-Cell Correlations','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');

%% Plot Histogram of all Quiet Cell-Cell Correlations by Genotype

WT_pcorrs = [];
scn1a_pcorrs = [];

for i = 1:length(data)
    
    lower_corr = [];
    lower_corr = tril(data(i).corrs_quiet, -1);
    lower_corr = reshape(lower_corr,[],1);
    lower_corr(lower_corr == 0) = [];
    
    if ismember(data(i).mouse(1),WT)
        WT_pcorrs = [WT_pcorrs; lower_corr];
    elseif ismember(data(i).mouse(1),scn1a)
        scn1a_pcorrs = [scn1a_pcorrs; lower_corr];
    end
    
end

% % Create x vectors for swarm plots
% WT_x = ones(length(WT_pcorrs),1);
% scn1a_x = 2*ones(length(scn1a_pcorrs),1);
% 
% % Plot all correlation coefficients as swarm charts with overlaid box plots
% swarmchart(WT_x,WT_pcorrs,'k');
% hold on
% swarmchart(scn1a_x,scn1a_pcorrs,'r');
% hold on
% g1 = repmat({'WT'},length(WT_x),1);
% g2 = repmat({'Scn1a'},length(scn1a_x),1);
% g = [g1; g2];
% boxplot([WT_pcorrs;scn1a_pcorrs],g,'Colors','w');
% title('Correlation Coefficients Between Ndnf+ Cells');
% xlabel('Genotype');
% ylabel('Correlation Coefficient');
% 
% WT
plot_WT_table = WT_pcorrs';

N = hist(plot_WT_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_WT_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_WT_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.5;

hold on;

% Scn1a

plot_scn1a_table = scn1a_pcorrs';

N = hist(plot_scn1a_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_scn1a_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_scn1a_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.5;

hold on;
plot([mean(plot_WT_table) mean(plot_WT_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(plot_scn1a_table) mean(plot_scn1a_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_WT = 'WT Mean Corr: ' + string(round(mean(plot_WT_table),2));
txt_scn1a = 'Scn1a+/- Mean Corr: ' + string(round(mean(plot_scn1a_table),2));
text((mean(plot_WT_table)-0.05),(max(N)+25),txt_WT,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_scn1a_table)-0.05),(max(N)-100),txt_scn1a,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cell Pairs','fontsize',14);
    title('Correlation Between Pairs of Ndnf-INs During Quiet Rest','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');

%% Plot Shuffled Histogram of all Cell-Cell Correlations by Genotype

WT_pcorrs_shuffle = [];
scn1a_pcorrs_shuffle = [];

for i = 1:length(data)
    
    lower_corr = [];
    lower_corr = tril(data(i).corrs_quiet_shuffle, -1);
    lower_corr = reshape(lower_corr,[],1);
    lower_corr(lower_corr == 0) = [];
    
    if ismember(data(i).mouse(1),WT)
        WT_pcorrs_shuffle = [WT_pcorrs_shuffle; lower_corr];
    elseif ismember(data(i).mouse(1),scn1a)
        scn1a_pcorrs_shuffle = [scn1a_pcorrs_shuffle; lower_corr];
    end
    
end

% % Create x vectors for swarm plots
% WT_x = ones(length(WT_pcorrs),1);
% scn1a_x = 2*ones(length(scn1a_pcorrs),1);
% 
% % Plot all correlation coefficients as swarm charts with overlaid box plots
% swarmchart(WT_x,WT_pcorrs,'k');
% hold on
% swarmchart(scn1a_x,scn1a_pcorrs,'r');
% hold on
% g1 = repmat({'WT'},length(WT_x),1);
% g2 = repmat({'Scn1a'},length(scn1a_x),1);
% g = [g1; g2];
% boxplot([WT_pcorrs;scn1a_pcorrs],g,'Colors','w');
% title('Correlation Coefficients Between Ndnf+ Cells');
% xlabel('Genotype');
% ylabel('Correlation Coefficient');
% 
% WT Neural Activity and Running
plot_WT_table = WT_pcorrs_shuffle';

N = hist(plot_WT_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_WT_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_WT_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.5;

hold on;

% Scn1a Neural Activity and Running

plot_scn1a_table = scn1a_pcorrs_shuffle';

N = hist(plot_scn1a_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_scn1a_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_scn1a_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.5;

hold on;
plot([mean(plot_WT_table) mean(plot_WT_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(plot_scn1a_table) mean(plot_scn1a_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_WT = 'WT Mean Corr: ' + string(round(mean(plot_WT_table),2));
txt_scn1a = 'Scn1a+/- Mean Corr: ' + string(round(mean(plot_scn1a_table),2));
text((mean(plot_WT_table)-0.05),(max(N)+25),txt_WT,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_scn1a_table)-0.05),(max(N)-150),txt_scn1a,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cell Pairs','fontsize',14);
    title('Shuffled Correlation Between Pairs of Ndnf-INs During Quiet Rest','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
%% Plot WT Real and Shuffled Cell-Cell Correlations during Quiet Rest

% WT Real Data
plot_real_table = WT_pcorrs';

N = hist(plot_real_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_real_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_real_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b1.FaceAlpha = 0.6;

hold on;

% Scn1a Neural Activity and Running

plot_shuffle_table = WT_pcorrs_shuffle';

N = hist(plot_shuffle_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_shuffle_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_shuffle_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',WT_color,'EdgeColor','none');
b2.FaceAlpha = 0.3;

hold on;
plot([mean(plot_real_table) mean(plot_real_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
plot([mean(plot_shuffle_table) mean(plot_shuffle_table)],[0 max(N)+15],'--','Color',WT_color,'LineWidth',2.5);
hold on;
txt_real = 'Real Mean Corr: ' + string(round(mean(plot_real_table),2));
txt_shuffle = 'Shuffled Mean Corr: ' + string(round(mean(plot_shuffle_table),2));
text((mean(plot_real_table)-0.05),(max(N)+15),txt_real,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_shuffle_table)-0.05),(max(N)-120),txt_shuffle,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cell-Cell Pairs','fontsize',14);
    title('WT Cell-Cell Correlations During Quiet Rest','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
%% Plot Scn1a Real and Shuffled Cell-Cell Correlations During Quiet Rest

% Scn1a Real Data
plot_real_table = scn1a_pcorrs';

N = hist(plot_real_table, [-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_real_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_real_table)*0.05*pdf(pd,x_g); % Plot pdf for data

figure;    
b1 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b1.FaceAlpha = 0.6;

hold on;

% Scn1a Shuffled Data

plot_shuffle_table = scn1a_pcorrs_shuffle';

N = hist(plot_shuffle_table,[-1:0.05:1]); % Create histogram of correlation coeffs
pd = fitdist(plot_shuffle_table.','normal'); % Fit probability dist to data
q_g = icdf(pd,[0.0013499 0.99865]); % Create inverse cumulative dist fxn from data
x_g = linspace(q_g(1),q_g(2)); % Create x-axis for histogram
y_g = numel(plot_shuffle_table)*0.05*pdf(pd,x_g); % Plot pdf for data

b2 = bar([-1:0.05:1],N,'FaceColor',scn1a_color,'EdgeColor','none');
b2.FaceAlpha = 0.3;

hold on;
plot([mean(plot_real_table) mean(plot_real_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
plot([mean(plot_shuffle_table) mean(plot_shuffle_table)],[0 max(N)+15],'--','Color',scn1a_color,'LineWidth',2.5);
hold on;
txt_real = 'Real Mean Corr: ' + string(round(mean(plot_real_table),2));
txt_shuffle = 'Shuffled Mean Corr: ' + string(round(mean(plot_shuffle_table),2));
text((mean(plot_real_table)-0.05),(max(N)+15),txt_real,...
    'FontSize',14,'HorizontalAlignment', 'right');
hold on;
text((mean(plot_shuffle_table)-0.05),(max(N)-140),txt_shuffle,...
    'FontSize',14,'HorizontalAlignment', 'right');

%hold on;plot(x_g,y_g,'g','LineWidth',2);
%     hold on;bar([-1:0.05:1],N_r,'FaceColor','w','EdgeColor','r')
%     hold on;plot(x_r,y_r,'r','LineWidth',2);
    xlabel('Correlation Coefficient','fontsize',14);
    ylabel('Number of Cell-Cell Pairs','fontsize',14);
    title('Scn1a Cell-Cell Correlations During Quiet Rest','fontsize',14);
    
ax = gca;
    
    set(gcf,'position',[100,100,540,400]);
    set(gca, 'FontName', 'Arial');
    
    
%% Calculate Mean Active Time Per Cell

WT_pct_act = [];
scn1a_pct_act = [];

for i = 1:length(data) % Loop over each recording
    
    for cell = 1:height(data(i).dF) % Loop over each cell
        
        %%% Define baseline fluoresence
        % Get local min of trace
        [min_v,min_i] = min(data(i).dF(cell,0.5*fs:end-0.5*fs));
        min_i = min_i + 0.5*fs;
        % Calculate average and SD for 500ms around min value
        mean_base = mean(data(i).dF(cell,min_i-0.5*fs:min_i+0.5*fs));
        sd_base = std(data(i).dF(cell,min_i-0.5*fs:min_i+0.5*fs));
        % Set threshold for transient as 3*SD of baseline
        t_thresh = 4*sd_base;

        %%% Calculate time spent over threshold
        active = data(i).dF(cell,:);
        active(active < t_thresh) = NaN;

        active_pct = sum(not(isnan(active)),2)/length(active); % pct of time active
        
        if ismember(data(i).mouse(1),WT)
            WT_pct_act = [WT_pct_act active_pct];
        elseif ismember(data(i).mouse(1),scn1a)
            scn1a_pct_act = [scn1a_pct_act active_pct];
        end
    
    end

end

% Create x vectors for swarm plots
WT_x = ones(length(WT_pct_act),1);
scn1a_x = 2*ones(length(scn1a_pct_act),1);

% Plot all correlation coefficients as swarm charts with overlaid box plots
swarmchart(WT_x,WT_pct_act,'MarkerFaceColor','k',...
    'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on
swarmchart(scn1a_x,scn1a_pct_act,'MarkerFaceColor',scn1a_color,...
     'MarkerFaceAlpha',0.4,'MarkerEdgeColor','none');
hold on
g1 = repmat({'WT'},length(WT_pct_act),1);
g2 = repmat({'Scn1a'},length(scn1a_pct_act),1);
g = [g1; g2];
bp = boxplot([WT_pct_act';scn1a_pct_act'],g,'Colors','k','BoxStyle','outline');
set(bp, 'LineWidth',1.5);
title('Percent of Time Active During Recording');
xlabel('Genotype');
ylabel('Percent time of recording');

set(gcf,'position',[100,100,540,400]);
set(gca, 'FontName', 'Arial');

%% Calculate Mean Active Time Per Cell During Rest

WT_pct_act = [];
scn1a_pct_act = [];

for i = 1:length(data) % Loop over each recording
    
    for cell = 1:height(data(i).dF) % Loop over each cell
        
        %%% Define baseline fluoresence
        % Get local min of trace
        [min_v,min_i] = min(data(i).dF(cell,0.5*fs:end-0.5*fs));
        min_i = min_i + 0.5*fs;
        % Calculate average and SD for 500ms around min value
        mean_base = mean(data(i).dF(cell,min_i-0.5*fs:min_i+0.5*fs));
        sd_base = std(data(i).dF(cell,min_i-0.5*fs:min_i+0.5*fs));
        % Set threshold for transient as 4*SD of baseline
        t_thresh = 4*sd_base;

        %%% Calculate time spent over threshold
        quiet_data = []; % array to store data during quiet periods
        for bout = 1:length(data(i).dF_quiet) % loop over each cell
            quiet_data = [quiet_data data(i).dF_quiet{bout}(cell,:)];
        end
        active = quiet_data;
        active(active < (mean_base + t_thresh)) = NaN;

        active_pct = sum(not(isnan(active)),2)/length(active) * 100; % pct of time active
        
        if ismember(data(i).mouse(1),WT)
            WT_pct_act = [WT_pct_act active_pct];
        elseif ismember(data(i).mouse(1),scn1a)
            scn1a_pct_act = [scn1a_pct_act active_pct];
        end
    
    end

end

% Create x vectors for swarm plots
WT_x = ones(length(WT_pct_act),1);
scn1a_x = 2*ones(length(scn1a_pct_act),1);

% Plot all correlation coefficients as swarm charts with overlaid box plots
swarmchart(WT_x,WT_pct_act,'MarkerFaceColor','k',...
    'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on
swarmchart(scn1a_x,scn1a_pct_act,'MarkerFaceColor',scn1a_color,...
     'MarkerFaceAlpha',0.4,'MarkerEdgeColor','none');
hold on
g1 = repmat({'WT'},length(WT_pct_act),1);
g2 = repmat({'Scn1a'},length(scn1a_pct_act),1);
g = [g1; g2];
bp = boxplot([WT_pct_act';scn1a_pct_act'],g,'Colors','k','BoxStyle','outline');
set(bp, 'LineWidth',1.5);
title('Percent of Time Active During Rest');
xlabel('Genotype');
ylabel('Percent time of recording');

set(gcf,'position',[100,100,540,400]);
set(gca, 'FontName', 'Arial');

%% Calculate Mean dF/F z-score for each Cell

WT_act = [];
scn1a_act = [];

for i = 1:length(data) % Loop over each recording
    
    quiet_data = [];
    for bout = 1:length(data(i).dF_quiet) % loop over each cell
        quiet_data = [quiet_data data(i).dF_quiet{bout}];
    end
    mean_F = mean(zscore(quiet_data,0,1),2);
    data(i).mean_zscore = mean_F;
    
    if ismember(data(i).mouse(1),WT)
        WT_act = [WT_act; mean_F];
    elseif ismember(data(i).mouse(1),scn1a)
        scn1a_act = [scn1a_act; mean_F];
    end

end

% Create x vectors for swarm plots
WT_x = ones(length(WT_act),1);
scn1a_x = 2*ones(length(scn1a_act),1);

% Plot all correlation coefficients as swarm charts with overlaid box plots
swarmchart(WT_x,WT_act,'MarkerFaceColor','k',...
    'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on
swarmchart(scn1a_x,scn1a_act,'MarkerFaceColor',scn1a_color,...
     'MarkerFaceAlpha',0.4,'MarkerEdgeColor','none');
hold on
g1 = repmat({'WT'},length(WT_act),1);
g2 = repmat({'Scn1a'},length(scn1a_act),1);
g = [g1; g2];
bp = boxplot([WT_act;scn1a_act],g,'Colors','k','BoxStyle','outline');
set(bp, 'LineWidth',1.5);
title('Mean dF/F During Rest');
xlabel('Genotype');
ylabel('Mean dF/F');

set(gcf,'position',[100,100,540,400]);
set(gca, 'FontName', 'Arial');