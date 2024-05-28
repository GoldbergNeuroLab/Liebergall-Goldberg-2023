%% Analysis of Intrinsic Properties from Current Steps (Tonic Spike)
% Script Author: Sophie Liebergall
% Updated: 1-26-23
% Purpose: Analysis of intrinsic properties from a typical current steps
% protocol (e.g. 20 pA current steps starting from -40 pA) calculating
% properties of single APs from first spike at tonic spiking (defined as
% first sweep in which >10 APs were fired) rather than rhe first spike
% during the rheobase sweep
% Inputs:
%   single .abf file with trace from classic current steps protocol in IC
% Outputs:
%   Vm (mV) 
%   Rm (mV) 
%   Time constant 
%   Rheobase (pA) 
%   AP Threshold (mV) 
%   AP Rise Time (ms) 
%   Max Rise Slope (mV/ms) 
%   AP Halfwidth (ms) 
%   AP Amplitude (mV) 
%   AHP Amplitude (mV) 
%   Sag (%) 
%   APs at rheobase 
%   Max instantaneous firing freq (Hz) 
%   Max steady state firing freq (Hz) 
% Dependencies:
%   - requires abfload.m and extrema.m to be in MATLAB search path

%% Import File
[data,si,h]=abfload(pathFileName);
% data = data frame configured as <data pts per sweep> by <number of chans> 
%   by <number of sweeps>.
% si = sampling interval in microseconds (us)
%       Note: sampling rate (Hz) = 1e6 / si (us)
% h = information on file (selected header parameters)
[~,fileName,~] = fileparts(pathFileName);%savepath = fullfile(saveRoot,saveDir1,saveDir2,fileName);

% Extract dimensions of the full data matrix
[i, num_chan, num_sweeps]=size(data); 
% i = number of samples in each sweep (i.e. recorded data points)
% num_chan = number of channels
% num_sweeps = number of sweeps in .abf file

% Specify channel of interest (assuming first channel is membrane voltage
% and second is current)
% User may need to change this value depending on channel of interest

if num_chan < 4
    v_chan = 1;
    i_chan = 2;
else    
    v_chan = 3; % Membrane voltage channel
    i_chan = 4; % Secondary output reading current channel
end
    
% Create vector with all indices in file
indices = 1:i;

% Create a vector with all indices converted to time points
time(indices,1) = (indices./(1e6./si));

% Set threshold for defining an action potential (in mV)
threshold = 0;

% Find on and off times for current steps
I_diff = diff(data(:,i_chan,1)); % extract current signal (channel 2) for first sweep
I_on = find(I_diff == min(I_diff),1); % get index when current goes on (minimum of dI/dt)
I_off = find(I_diff == max(I_diff),1); % get index when current goes off (maximum of dI/dt)
% Note: this assumes that this time is the same for each sweep and that the
% first sweep injects negative current

%% Create Table of Spikes for Each Sweep

spike_table = zeros(num_sweeps,4); % create empty table to store spikes for each sweep
% Formatted as [Sweep] [Absolute Current (pA)] [Difference from Holding (pA)] [Number of Spikes] 


spike_i = {}; % create empty cell array to store indices of spikes for each sweep (to calc ISIs)

for sweep = 1:num_sweeps % loop through each sweep in the .abf file
    v = data(I_on:I_off+100,v_chan,sweep); % extract voltage data (channel 1) for current sweep
    %v = smooth(v, 1); % smooth data to avoid double counting APs (seemed
    %to lead to undercounting of APs)
    
    k = find(v > threshold);
    if isempty(k)
        spike_is = [];
    else
        spike_is = [k(1)];
        for w = 2:length(k)
            if k(w) ~= k(w-1)+1
                spike_is = [spike_is, k(w)];
            end
        end
    end    
    spike_i{sweep} = spike_is; % get indices of local maxima above threshold (i.e. spikes)
    
    I = data(int16((I_on+I_off)/2),i_chan,sweep); % get current value for sweep
    I_hold = data(int16(I_on-100),i_chan,sweep);
    
    spike_table(sweep,1) = sweep; % record sweep #
    spike_table(sweep,2) = I; % record absolute current for sweep
    spike_table(sweep,3) = I - I_hold; % record current difference from hold for sweep
    spike_table(sweep,4) = size(spike_i{sweep},2); % record number of spikes for a sweep
    
end;

%% Max Steady State Firing Frequency

% Get the sweep with the greatest number of spikes from the spike table
[max_spikes,max_sweep] = max(transpose(spike_table(:,4))); 

% Calculate max steady state firing frequency
max_SSFF = max_spikes./((I_off - I_on)./1e6.*si); % in Hz

%% Max Instantaneous Firing Frequency

% Create cell array with inter-spike intervals (ISIs)
ISIs = {};

for sweep = 1:num_sweeps % loop through each sweep
    ISI_sweep = []; % create vector to store ISIs for each sweep
    if spike_table(sweep,4) > 1
        for int = 2:spike_table(sweep,4); % loop through each ISI
            % Calculate ISI and append to ISI_sweep vector
            ISI_sweep = [ISI_sweep (spike_i{sweep}(int) - spike_i{sweep}(int-1))];
        end
    ISIs{sweep} = ISI_sweep; % Add ISIs for a sweep to ISIs cell array
    end
end

sorted_ISIs = sort(cell2mat(ISIs)); 

% Use while loop to ensure that Max IFF is within reasonable range (<2000
% Hz) - highest recorded max IFF for PV-INs around 1700 Hz
n = 1;
max_IFF = 2001;
while max_IFF > 2000;
    max_IFF = 1./(sorted_ISIs(n)./1e6.*si); % calcualte max IFF from min ISI (Hz)
    n = n + 1;
end
    
%% Rheobase (and spikes at rheobase)

% Get sweep number of rheobase (rheo_sweep)
spike_sweeps = spike_table(spike_table(:,4) > 0);
rheo_sweep = spike_sweeps(1);
% Get current injection (pA) during rheo_sweep sweep
rheobase = spike_table(rheo_sweep,3);

spikes_at_rheobase = spike_table(rheo_sweep,4);

%% FF at 2x Rheobase

% Get sweep number where current injection = 2*rheobase
double_rheo_sweep = spike_table(spike_table(:,3) >= rheobase*2);
if ~isempty(double_rheo_sweep)
    double_rheo_sweep = double_rheo_sweep(1);
    % Get firing frequency during 2*rheobase sweep
    double_rheo_FF = spike_table(double_rheo_sweep,4)./((I_off - I_on)./1e6.*si);
else
    double_rheo_FF = max_SSFF;
end

%% Tonic Sweep

if max(spike_table(:,4)) < 10
    tonic_sweep = rheo_sweep;
else
    ton_sweeps = spike_table(spike_table(:,4) > 9);
    tonic_sweep = ton_sweeps(1);
end

%% Resting membrane potential

% Get sweep with 0 current injection (i.e. first non-negative current
% injection)
non_neg_sweeps = spike_table(spike_table(:,2) >= 0);
zero_sweep = non_neg_sweeps(1);

% Calculate mean voltage when injecting zero current
Vm = mean(data(I_on:I_off,v_chan,zero_sweep));

%% Input Resistance
% Calculated from first trace (assumed to be negative I injection)

% Calcualte mean voltage before current step
V1 = mean(data(1:I_on,v_chan,1));
% Calclate mean voltage in second half of current step
V2 = mean(data(((I_on+I_off)./2):I_off,v_chan,1));

% Calculate mean current before current step
I1 = mean(data(1:I_on,i_chan,1));
% Calclate mean current in second half of current step
I2 = mean(data(((I_on+I_off)./2):I_off,i_chan,1));

Rm = (V2 - V1)./(I2 - I1).*1000; % Calculate Rm in mOhms

%% Membrane time constant
% Calcualte expontential fit for 20 ms from 1st sweep (assuming negative
% current injection)

% Set fit duration in ms
fit_dur = 20;

% Create time bin for calculating membrane time constant from first sweep
fit_data = data(I_on:(I_on+(fit_dur.*1e3./si)),v_chan,1);
fit_data_time = transpose((1:fit_dur*1e3./si+1)./1e3.*si);

% Specify fit type as single expontential with c term
ft = fittype('a*exp(-b*t) + c','indep','t'); 
start_pt = [6,0.1,-65]; % Set start point for fit coefficients (from literature)

f = fit(fit_data_time,fit_data,ft,'StartPoint',start_pt); % Fit curve
% figure;
% plot(f,fit_data_time,fit_data)
membrane_time_constant = 1/f.b;

%% Sag Percentage
% % Sag (produced by Ih) is a ratio of the steady state relative to the 
% maximal hyperpolarization produced in response to a negative current 
% injection from resting membrane potential. 
% Note: calcualted from sweep 1 which is presumed to be a negative current
% injection

sag_min = min(data(I_on:I_off,v_chan,1)); % Get min voltage during negative current inj.

% Calculate steady state as voltage 3/4 of way through current step
steady_index = int16((I_off - I_on)*0.75 + I_on);
sag_steady = data(steady_index,v_chan,1);

sag = (1 - sag_steady/sag_min)*100;

%% AP Threshold
% Calculated via 2 methods

% Isolate voltage data from current step on to peak of 1st AP at rheobase
% AP1_v is onset of current step to 5 ms after peak of first AP

endtime = round(0.003*(1e6./si));
if length(spike_i{tonic_sweep}) == 1
    otherside = 0.1*(1e6./si);
else
    otherside = spike_i{tonic_sweep}(2) - spike_i{tonic_sweep}(1) - endtime;
end

AP1_v = data(I_on:(I_on+spike_i{tonic_sweep}(1)+otherside),v_chan,tonic_sweep);

AP1_v_d = gradient(AP1_v)./(si*1e-3); % get 1st derivative of voltage data (mV/ms)
AP1_v_d2 = gradient(AP1_v_d)./(si*1e-3); % get 2nd derivative of voltage data (mV/ms)

%%% Method 1: Where 1st derivate = 20 mV/ms
AP_thresh1i = find(AP1_v_d(100:end) >= 20);
if isempty(AP_thresh1i)
    AP_thresh1i = 0;
else
    AP_thresh1i = AP_thresh1i(1)+100; % Get index where dV/dt >/ 20 mV/ms
    AP_thresh1 = AP1_v(AP_thresh1i);
end

%%% Method 2: Where d2V/dt is 10% of max value
max_d2_10pct = max(AP1_v_d2)*0.1; % get 10th percentile of max d2V/dt value
AP_thresh2i = find(AP1_v_d2(100:end) >= max_d2_10pct);
if isempty(AP_thresh2i)
    AP_thresh2i = 0;
else
    AP_thresh2i = AP_thresh2i(1)+100;
    AP_thresh2 = AP1_v(AP_thresh2i);
end

%% AP Peak

[AP_peak,AP_peaki] = max(AP1_v);

%% AP Rise Time
% Calculated as difference between time at AP peak and time at AP threshold
AP_rise_time = (AP_peaki - AP_thresh1i)./1e3.*si; %in ms

%% Max AP Rise Slope
% Calculated as max of 1st derivative between AP threshold and peak in
% mV/ms
max_rise_slope = max(AP1_v_d(AP_thresh1i:AP_peaki));

%% AP Amplitude

AP_amplitude = AP_peak - AP_thresh1;

%% AP Halfwidth
% AP half-width (AP ½-width) is defined as the width of the AP (in ms) at 
% half-maximal amplitude, calculated using AP threshold and AP peak. 

% Find index of 1/2 AP amplitude on upstroke
AP_half_up = find(AP1_v >= (AP_amplitude*0.5+AP_thresh1));
AP_half_up = AP_half_up(1);

% Find index of 1/2 AP amplitude on downstroke
AP_half_down = find(AP1_v(AP_peaki:end) <= (AP_amplitude*0.5+AP_thresh1));
AP_half_down = AP_half_down(1) + AP_peaki - 1;

% Calculate AP half_width
AP_half_width = (AP_half_down - AP_half_up)./1e3*si;

%% After Hyperpolarization (AHP) Amplitude
% Calculated as difference between AP threshold and AHP trough V (in mV)

[AHP_min,AHP_min_time] = min(AP1_v(AP_peaki:end));
AHP_amplitude = AP_thresh1 - AHP_min;

%% Latency to 1st Spike

% Isolate voltage data from current step on to peak of 1st AP at rheobase
% AP1_v is onset of current step to 5 ms after peak of first AP

if length(spike_i{rheo_sweep}) == 1
    otherside = 0.1*(1e6./si);
else
    otherside = spike_i{rheo_sweep}(2) - spike_i{rheo_sweep}(1) - endtime;
end

AP1_vr = data(I_on:(I_on+spike_i{rheo_sweep}(1)+otherside),v_chan,rheo_sweep);

AP1_vr_d = gradient(AP1_vr)./(si*1e-3); % get 1st derivative of voltage data (mV/ms)

%%% Method 1: Where 1st derivate = 20 mV/ms
APr_thresh1i = find(AP1_vr_d(100:end) >= 20);
if isempty(APr_thresh1i)
    APr_thresh1i = 0;
else
    APr_thresh1i = APr_thresh1i(1)+100; % Get index where dV/dt >/ 20 mV/ms
    APr_thresh1 = AP1_vr(APr_thresh1i);
end

latency_to_AP1 = APr_thresh1i./1e3*si;

%% Format Intrinsic Properties as Table
% Create MATLAB table with intrinsic properties

IntrinsicProperty = ["Vm (mV)";"Rm (mOhms)"; "Time Constant";"Rheobase (pA)"; ...
    "AP Threshold 1 (mV)";"AP 2(mV)";"AP Rise Time (ms)";"Max Rise Slope (mV/ms)"; ...
    "AP Halfwidth (ms)";"AP Peak (mV)";"AP Amplitude (mV)"; "AHP Amplitude (mV)";...
    "Sag (%)";"APs at Rheobase";"Max IFF (Hz)";"Max SSFF (HZ)";"2X Rheo FF";"Latency to 1st Spike (ms)"];

Values = [Vm; Rm; membrane_time_constant; rheobase; ...
    AP_thresh1; AP_thresh2; AP_rise_time; max_rise_slope; ...
    AP_half_width; AP_peak; AP_amplitude; AHP_amplitude;...
    sag; spikes_at_rheobase; max_IFF; max_SSFF; double_rheo_FF;...
    latency_to_AP1];

intrinsic_properties = table(IntrinsicProperty, Values);

%% Create Phase Plot of Rheobase Action Potential
% 
 x_mV = AP1_v(1:(length(AP1_v)-endtime));
 y_dVdt = gradient(x_mV)./(si*1e-3);
% 
fig = figure;
plot(x_mV,y_dVdt)
xlim([-80 40])
ylim([-150 350])
xlabel('Membrane Potential (mV)','FontName','Arial');
ylabel('Membrane Potential Change (mV/ms)','FontName','Arial');
title(fileName,'FontName','Arial');
% 
% % Save the phase plot to analysis folder
 saveRoot = '\\ressmb01.research.chop.edu\goldberg_lab\Sophie Liebergall\analysis\electrophysiology\ndnf_scn1a_intrinsic_properties\ndnf_scn1a_P35-56\';
 saveDir1 = 'WT\';
 saveDir2 = 'phase_plots\';
 savepath = fullfile(saveRoot,saveDir1,saveDir2,fileName);
% 
 saveas(fig,savepath,'pdf');
 saveas(fig,savepath,'png');
 
 close

%% Create IF Plot

%plot(spike_table(rheo_sweep:end,2),spike_table(rheo_sweep:end,3));

% Save the spike table to specified file path so can be accessed to create
% average IF curves
saveDir2 = 'IF_curves\';

%writematrix(spike_table,savepath,'FileType','text');

%% Create Plot of 3rd Derivative of AP

AP1_v_d3 = gradient(AP1_v_d2)./(si*1e-3); % get 3rd derivative of voltage data (mV/ms)
x_t = (1:length(AP1_v_d3))./(1e3./si);

% Plot 1ms before and 3ms after
xlim1 = AP_thresh1i./(1e3./si) - 1;
xlim2 = xlim1 + 3;

fig2 = figure;
plot(x_t, AP1_v_d3);
xlim([xlim1 xlim2])
%ylim([-150 350])
xlabel('Time (ms)','FontName','Arial');
ylabel('d3V/dt3 (mV/ms^3)','FontName','Arial');
title(fileName,'FontName','Arial');

% Save a plot of the 3rd derivative of the AP to observe a stationary
% inflection
saveDir3 = 'AP_3deriv\';
savepath = fullfile(saveRoot,saveDir1,saveDir3,fileName);

 saveas(fig2,savepath,'pdf');
 saveas(fig2,savepath,'png');
 
 close
