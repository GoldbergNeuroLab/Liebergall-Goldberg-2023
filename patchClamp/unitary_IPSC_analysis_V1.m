%% Analysis of Unitary IPSCs (Outward)
% Script Author: Sophie Liebergall
% Updated: 1-13-2023
% Purpose: Analysis of unitary IPSCs
% Inputs:
%   single .abf file with VC trace of a unitary IPSC
% Outputs:
%   Amplitude, synaptic latency, 10-90 rise time, 10-90 decay time, 
% Dependencies:
%   - requires abfload.m and extrema.m to be in MATLAB search path

%% Import File
[data,si,h]=abfload(pathFileName);
% data = data frame configured as <data pts per sweep> by <number of chans> 
%   by <number of sweeps>.
% si = sampling interval in microseconds (us)
%       Note: sampling rate (Hz) = 1e6 / si (us)
% h = information on file (selected header parameters)

% Extract dimensions of the full data matrix
[i, num_chan, num_sweeps]=size(data); 
% i = number of samples in each sweep (i.e. recorded data points)
% num_chan = number of channels
% num_sweeps = number of sweeps in .abf file

% Specify channel of interest (assuming first channel is membrane voltage
% and second is current)
% User may need to change this value depending on channel of interest

v_chan_stim = 1;
i_chan_stim = 3;
i_chan_IPSC = 2;
v_chan_IPSC = 4;
    
% Create vector with all indices in file
indices = 1:i;

% Create a vector with all indices converted to time points
time(indices,1) = (indices./(1e6./si));

%% Get Index at Peak of AP

% Create table to store indiced of AP peak for each sweep
peak_table = zeros(num_sweeps,3);

for sweep = 1:num_sweeps
    [AP_max,I_AP] = max(data(:,v_chan_stim,sweep));
    peak_table(sweep, 1) = sweep;
    peak_table(sweep, 2) = AP_max;
    peak_table(sweep, 3) = I_AP;   
end

%% Create average trace from aligned IPSCs

aligned_traces = {};

baseline_win = 0.1*(1e6./si);
% added bc of inward current from electrical connection during presyn AP
upstroke_time = 0.001*(1e6./si); 

for sweep = 1:num_sweeps
    if peak_table(sweep,2) > 0
        AP_peaki = peak_table(sweep,3);
        starti = AP_peaki - baseline_win - upstroke_time;
        IPSC = data(starti:AP_peaki+(.2*(1e6./si)),i_chan_IPSC,sweep);
        
        aligned_traces{sweep} = IPSC;
    end
end

% Generate average trace
aligned_tracesm = cat(3,aligned_traces{:});
average_trace = mean(aligned_tracesm,3);
average_trace_s = smoothdata(average_trace);

%% Calculate average amplitude

[max_val, amp_I] = max(average_trace);
baseline = mean(average_trace(1:(.1*(1e6./si))));
amplitude = max_val - baseline;

%% Calculate time to peak

time_to_peak = (amp_I - baseline_win - upstroke_time)/(1e3./si); % in ms

%% Calculate 10-90 rise time

rise10 = find(average_trace((baseline_win + upstroke_time):end) >= (baseline+(amplitude*.1)),1);
rise90 = find(average_trace((baseline_win + upstroke_time):end) >= (baseline+(amplitude*.9)),1);

risetime10_90 = (rise90-rise10)/(1e3./si); % in ms

%% Calculate 10-90 decay time

decay10 = find(average_trace(amp_I:end) <= (average_trace(amp_I)-amplitude*.1),1);
decay90 = find(average_trace(amp_I:end) <= (average_trace(amp_I)-amplitude*.9),1);

decaytime90_10 = (decay90-decay10)/(1e3./si); % in ms

%% Calculate SD & RMS noise of baseline signal

rms_val = rms(average_trace(1:baseline_win)) - baseline;
sd = std(average_trace(1:baseline_win));

%% Format IPSC Properties as Table
% Create MATLAB table with intrinsic properties

IPSCProperties = ["Amplitude (pA)";"Time to Peak (ms)";"Rise Time (ms)"; ...
    "Decay Time (ms)";"RMS";"SD"];

Values = [amplitude; time_to_peak; risetime10_90; decaytime90_10; rms_val; sd];

IPSC_properties = table(IPSCProperties, Values);

