%%% Align Pupillometry, Locomotion, and 2P Imaging Neural Data

% Author: Sophie Liebergall (adapted from Ala Somarowthu)
% Updated: 4/10/23
% Purpose: Align behavioral and neural 2P data
% Inputs:
% - suite2p output folder
% - .xlsx file with locomotion data from Neurotar software
% - .csv file with pupillometry data from DeepLabCut

clc;
clear all;
close all;

cur_dir=pwd;

data_dir='D:\two_photon\processed\ndnf_cre_arousal';

save_dir='D:\two_photon\processed\ndnf_cre_arousal';

dir_list=dir([data_dir '\*FOV*']);

for d_ind=1:length(dir_list)%1:length(dir_list)
    d_ind
    close all;
    clearvars -except fs cur_dir dir_name dir_list d_ind data_dir save_dir
    dir_n=[data_dir '\' dir_list(d_ind).name]; % path to file
    save_dir_n = fullfile(save_dir,dir_list(d_ind).name);
    mkdir(save_dir_n);
    
    %%% Loading suite2p file
    load (fullfile(dir_n,'processed_2p\suite2p\plane0\Fall.mat'));
    ind_slash=strfind(dir_n,'\');
    
    %%% Getting fs from xmlfile

    fn=dir([dir_n,'\processed_2P\*.xml']);

    Dm=readstruct(fullfile(dir_n,'\processed_2P\',fn(1).name));

    for i=1:length(Dm.Sequence.Frame(1).PVStateShard.PVStateValue)

        if Dm.Sequence.Frame(1).PVStateShard.PVStateValue(i).keyAttribute=="framePeriod"

            fs=1/Dm.Sequence.Frame(1).PVStateShard.PVStateValue(i).valueAttribute;

        end

    end
    
    %%% Read in pupil data
    temp_p=dir([dir_n, '\pupil\*.csv']);
    %temp_p=dir([dir_n, '\pupil\*.xlsx' dir_n(ind_slash(end)+1:ind_slash(end)+12) upper(dir_n(ind_slash(end)+13:ind_slash(end)+16)) '*.mat']);
    
    fname=fullfile([dir_n, '\pupil\'],temp_p(1).name);
    pupil_data = csvread(fname,4,1); % Start at 4th row
    
    % x1,y1 = north
    % x2,y2 = east
    % x3, y3 = south
    % x4, y4 = west
    pupil_vars = {'x1','y1','likelihood1','x2','y2','likelihood2',...
        'x3','y3','likelihood3','x4','y4','likelihood4'};
    
    pupil_table = array2table(pupil_data,'VariableNames',pupil_vars);
    
    % Calculate pupil diameter
    north = [pupil_table.x1, pupil_table.y1];
    east = [pupil_table.x2, pupil_table.y2];
    south = [pupil_table.x3, pupil_table.y3];
    west = [pupil_table.x4, pupil_table.y4];
    
    % Calculate euclidian distance between N/S and E/W points
    diam1 = sqrt((south(:,1)-north(:,1)).^2 + (south(:,2)-north(:,2)).^2);
    diam2 = sqrt((west(:,1)-east(:,1)).^2 + (west(:,2)-east(:,2)).^2);
    avg_diam = mean([diam1,diam2],2); % Find average of two diameters
    
    diam_masked = avg_diam; % Remove low condfidence values
    diam_masked(pupil_table.likelihood1 < 0.9) = NaN;
    
    pupil_table.diameter = diam_masked; % Add as column to pupil table
    
    pupil_fs = 100; % Set sample rate for pupil data in Hz
    
    %%% Loading running data
    % Add running bout scoring parameters
    speed_th=10;
    baseline_win=5;
    win1=4;
    win2=4;
    win_l=1;
    win1_m=15;
    win2_m=5;
    win_m=5;
    
    % Get path to running data
    fname=fullfile([dir_n, '\running\'] ,'track.xlsx');
    
    % Reading TTL pulse for start of imaging scan
    [num,str]=xlsread(fname,'Pp_Data');
    TTL=str2double(str(:,end-2));
    TTL_running=double(TTL==1);
    %TTL_whisk_puff=double((TTL==1000));
    
    % Aligning resampling data
    [~,st] = pulsewidth(TTL_running); % Gets start of 2P imaging scan
    win=[round(st(1)):length(num)]; % Window during imaging
    
    % Get MHC Tracker time stamps
    time_stamp=num(:,3);
    time_stamp=time_stamp(win)-time_stamp(round(st(1)));
    
    % Get instantenous running speed from MHC tracker
    inst_speed=num(:,12); 
    inst_speed=inst_speed(win);
    
    % Resample running & pupil data to match sampling rate of 2P imaging (30 Hz)
    [inst_speed_r,tr]=resample(inst_speed,time_stamp,fs);
    pupil_time_stamp = (1:height(pupil_table))./pupil_fs;
    [pupil_dia_f]=resample(diam_masked,pupil_time_stamp,fs);

    % Trim running data to length of 2P recording
    New_d=min([length(inst_speed_r),size(F,2),length(pupil_dia_f)]);

    inst_speed_r(New_d+1:end)=[];
    pupil_dia_f(New_d+1:end)=[];
    tr(New_d+1:end)=[];
    F(:,New_d+1:end)=[];
    Fneu(:,New_d+1:end)=[];
    spks(:,New_d+1:end)=[];  
    
    % Computing dF_F0
    [dF_F0_g]= dipoppa_neuropil_sub_1ch_edit(F,Fneu,find(iscell(:,1)==1),fs);

    [bb,aa]=butter(10,[2/fs]);%4HZ before
    inst_speed_r=filtfilt(bb,aa,inst_speed_r); % smooth running data
    pupil_dia_noblinks = filloutliers(pupil_dia_f,"next"); % remove blinks from pupil data (outliers)
    pupil_dia_f_r=medfilt1(pupil_dia_noblinks,10); % smooth with median filter    

    cd (save_dir_n)
    save Data.mat dF_F0_g Fneu fs inst_speed_r ops pupil_dia_f_r tr
    % % save Data.mat dF_F0_r dF_F0_g M_r M_g inst_speed_r dF_F0 inst_trace dF_F0_r_trace dF_F0_g_trace
    cd (cur_dir)
    
end
    
