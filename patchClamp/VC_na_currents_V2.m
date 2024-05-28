
%% Na+ current activation and fast inactivation from nucleated macropatch recordings

% Author: Sophie Liebergall
% Updated: 5/10/2023
% Adapted from Julie Merchant's code for VC recordings in heterologous
% systems

% Inputs: 
    % single .abf file of Na+ current activation protocol in voltage clamp
        % typically 300 ms depolarizing voltage steps from -80 mV to +55 mV in 5 mV increments

% User input per file:
    % fig 1 - click time bounds within which to identify peak current (exclude things related to compensation)
    % cell capacitance (TBD: read table with capacitance values calculated from separate script)
    % fig 4 - max voltage step to include in the fit (aka when it reaches  max conductance)

% Outputs:
    %Two results tables in workspace for the .abf file run
        %Activation_results has capacitance (input), V_hold, I_hold, peak absolute current (w/ leak subtraction), Vhalf activation, slope factor k activation)
        %Activation_results_per_step has 7 columns:
            % V_steps (column1)
            % Peak currents for IV (not leak subtracted) (column2)
            % Peak current density for normalized IV (not leak subtracted) (column3)
            % Time to peak (column4)
            % Normalized conductance (G/Gmax), set to 1 once it reaches max (column5)
            % Tau fast adn slow inactivation (probably only accurate ~peak current activation sweeps) (columns6&7)

% Dependencies:
    % abfload_HEK.m in MATLAB search path 
    % .abf file(s) to be analyzed in MATLAB search path 

%% Import File

close all
clearvars 

% Set indices for current and voltage channels
I_chan = 3;
V_chan = 2;

% Get the path of .abf file of interest and folder it resides in 
[fileName,path] = uigetfile('*.abf', 'abf files');
pathFileName = [path fileName];

% Make this the default path
cd(path)

% Load the data 
[data,si,h]=abfload(pathFileName);
    % data = data frame configured as <data pts per sweep> by <number of channels> by <number of sweeps>.
    % si = sampling interval in microseconds (us)
    %       Note: sampling rate (Hz) = 1e6 / si (us)
    % h = information on file (selected header parameters)

% Extract dimensions of full data matrix / set up loaded data 
[data_points, num_chan, num_sweeps]=size(data);
    % data_points = number of samples in each sweep (i.e. recorded data points) %i in Sophie's 
    % num_chan = number of channels
    % num_sweeps = number of sweeps in .abf file

% Create a vector with all indices in file and convert to time points
    indices = 1:data_points;
    time(indices,1)=(indices./(1e6./si)); % divide datapoints by sampling rate (Hz) to give time in s for each datapoint
    time_ms = time.*1000;   % gives time in ms for each datapoint

% Create 2D matrices from full data matrix
    matrix_V = zeros(data_points,num_sweeps);
    matrix_I = matrix_V;

    for i = 1:num_sweeps
        matrix_V(:,i) = data(:,V_chan,i);        %creates matrix of VOLTAGE for each sweep (column)
        matrix_I(:,i) = data(:,I_chan,i);        %creates matrix of CURRENT for each sweep (column)
    end

%% Identify holding voltage and current; voltage steps

% Identify holding voltage and current
V_hold = mean(matrix_V(1:100,1)); %find avg voltage in mV w/in time window in 1st sweep (hardcoded before pulse, b/w 2.1-4.2 ms)
I_hold = mean(mean(matrix_I(1:100,:))); %find avg current in pA w/in time window (hardcoded same as above) for each sweep, then avg all of those 

% Identify when voltage step goes on/off
    % assumes pulse time is the same for each sweep, and that 1st sweep is depolarizing from V_hold
V_diff = diff(matrix_V); %extract dV/dt
V_on = find(V_diff == max(V_diff),1); %get index when voltage step goes on
V_on_t = time_ms(V_on,1); %time step goes on in ms
V_off = find(V_diff == min(V_diff),1); %get index when voltage step goes off
V_off_t = time_ms(V_off,1); %time step goes off in ms 

% Loop to identify voltage step for each sweep
V_steps = zeros(1,num_sweeps); %preallocate
for i = 1:num_sweeps
    boolean = matrix_V(:,i)>mean(matrix_V(:,i)); %for each sweep, makes logical array for whether each value is > mean of whole sweep (essentially putting 1s when step is on)
    V_steps(i) = round(mean(matrix_V(boolean,i))); %for each sweep, find mean of all datapoints corresponding to when step is on) 
end

%% Plot and identify peak transient currents and time to peak 

% Plot transient Na+ currents
figure
plot(indices,matrix_I); %plots current traces vs datapoints
    xlim([0 1500]); %plots b/w datapoints translating to b/w ~7-15ms
    xticks(250:50:500)
    xticklabels({'0','1.5','3','4.5','6','7.5'});
    title('Transient Na+ currents');
    xlabel('Time (ms)');
    ylabel('current (pA)');
    %fig1name = ['Fig1 Currents ' num2str(fileName(1:end-4)) '.png'];
    %saveas(gcf,fig1name);
    
% Define x axis range w/in which to find peak [INSTRUCTIONS: click twice on plot to define bounds, then double click outside of plot to return to script]
[xi,yi] = getpts;
range_start = round(xi(1));
range_end = round(xi(2));

% Preallocate 
peak_I = zeros(1,num_sweeps);
range_peak_index = peak_I;
peak_index = peak_I;
peak_I_leaksub = peak_I;
peak_index_time = peak_I;
time_to_peak = peak_I;

% Loop to find peak transient currents (w/ and w/o leak subtraction) and time to peak 
for i = 1:num_sweeps       
    [peak_I(i), range_peak_index(i)] = min(matrix_I(range_start:range_end,i));    %find peak current (min) per sweep in time range and index of peak w/in that time range
    peak_index(i) = range_start + range_peak_index(i);     %find datapoint (overall index) corresponding to each peak per sweep
    peak_I_leaksub(i) = matrix_I(peak_index(i),i)-I_hold;  %subtract I_hold from current value corresponding to each peak per sweep
    peak_index_time(i) = time_ms(peak_index(i));  %find time in ms of each peak
    time_to_peak(i) = peak_index_time(i)-(V_on_t);   %find time in ms b/w V_on_t and time of peak per sweep
end

peak_I_abs = min(peak_I_leaksub);   %leak-subtracted max current size --> use this for exclusion later if >7 or 10 nA

%% Plot I-V; find and plot current density

% Plot I-V curve
figure
    plot(V_steps, peak_I,'k-o','LineWidth',1.5,'MarkerFaceColor','w');  
        title('I-V curve')
        xlim([-80 60]);
        xticks(-80:10:60);
        xlabel('voltage (mV)');
        ylabel('current (pA)');
        %fig2name = ['Fig2 IV ' num2str(fileName(1:end-4)) ' .png'];
        %saveas(gcf,fig2name);
        
% Get cell capacitance (proxy for cell size) from user input 
prompt = {'Enter cell capacitance in pF'}; %[INSTRUCTIONS: enter value determined elsewhere - Clampfit or VC_HEK_capacitance script]
answer = inputdlg(prompt);
capacitance = str2double(answer{1});

% Find current density
peak_I_density = peak_I/capacitance;  %normalize peak current size to cell capacitance

% Plot current density vs voltage
figure
    plot(V_steps,peak_I_density,'b-o','LineWidth',1.5,'MarkerFaceColor','w'); 
        xlim([-60 40]);
        xticks(-60:10:40);
        title('current density');
        xlabel('voltage (mV)');
        ylabel('current density (pA/pF)');
        %fig3name = ['Fig3 IV density ' num2str(fileName(1:end-4)) ' .png'];
        %saveas(gcf,fig3name);

%% Calculate and plot conductance, and normalize to Gmax and plot        

% Loop to find conductance per voltage step
G = zeros(1,num_sweeps);
Erev = 73.45166;  %hardcoded based on Erev calculated from Cl-F internal and standard external solutions
for i = 1:num_sweeps
    G(i) = peak_I(i)/(V_steps(i)-Erev);  %finds conductance using V=IR, G=1/R, G=I/deltaV
end

% Plot conductance vs voltage
figure
    plot(V_steps,G,'-o','LineWidth',1.5);
        title('conductance');
        xlim([-80 60]);
        xticks(-80:5:60);
        xlabel('voltage (mV)');
        ylabel('conductance (G)');
        grid on 
        %fig4name = ['Fig4 G ' num2str(fileName(1:end-4)) ' .png'];
        %saveas(gcf,fig4name);
        
% Find Gmax based on max voltage to include in the fit 
prompt = {'Enter max voltage (mV) to include in the fit'};
answer = inputdlg(prompt);
max_V = str2double(answer{1});
max_V_sweep = sum(V_steps<max_V+2); %find sweep of max voltage to include in the fit
Gmax = max(G(1:max_V_sweep));
Gmax_index = find(G == Gmax);  %find sweep of Gmax

% Loop to normalize G to Gmax for each sweep
norm_G = ones(1,num_sweeps); %preallocate WITH ONES to over-write in norm_G after reaching Gmax 
for i = 1:max_V_sweep 
    norm_G(:,i) = G(i)/Gmax;
end

% Plot normalized G
figure
    plot(V_steps,norm_G,'c-o','LineWidth',1.5);
        title('voltage dependence of activation');
        xlim([-80 40]);
        xticks(-80:10:40);
        xlabel('voltage (mV)');
        ylabel('normalized conductance (G/Gmax)');
        %fig5name = ['Fig5 normalized G ' num2str(fileName(1:end-4)) ' .png'];
        %saveas(gcf,fig5name);

%% Fit Boltzmann function to data and plot

% Perform Boltzmann fit
xval = transpose(V_steps);
yval = transpose(norm_G);
sig = fittype('1/(1+(exp((a-x)./b)))');
[boltz_fit,nnn] = fit(xval,yval,sig,'StartPoint',[1,1]);
Vhalf = boltz_fit.a;
k = boltz_fit.b;

% Plot normalized conductance datapoints w/ Boltzmann fit and Vhalf 
figure
    plot(boltz_fit,'k');
    hold on 
    plot(xval,yval,'.k','MarkerSize',15);
    hold on
    plot(Vhalf,0.5,'mp','MarkerSize',12,'MarkerFaceColor','m');
        legend('Boltzmann fit','G/Gmax','Vhalf','Location','northwest');
        title('voltage dependence of activation')
        xlim([-80 40]);
        xticks(-80:10:40);
        xlabel('voltage (mV)');
        ylabel('normalized conductance (G/Gmax)');
        %fig6name = ['Fig6 activation ' num2str(fileName(1:end-4)) ' .png'];
        %saveas(gcf,fig6name);
        
%% Find tau fast inactivation

% Preallocate matrices
tau_fast = zeros(1,num_sweeps);
tau_slow = zeros(1,num_sweeps);

% Loop to perform double exponential fit on inactivation decay per sweep
ft = fittype('exp2');
opts = fitoptions('Method','NonlinearLeastSquares');
    for i = 1:num_sweeps
        tau_startx_index = peak_index(i)+6;     %starts fit 6 datapoints (~0.18 ms) after peak for each sweep
        tau_endx_index = 1000;   %ends fit at 30 ms for each sweep (could alter this and see if you get better fits)
        t = time_ms(tau_startx_index:tau_endx_index);   %must input time in ms in order to get output tau in ms 
        j = matrix_I(tau_startx_index:tau_endx_index,i)-matrix_I(1000,i);
        [fitresult, gof] = fit(t, j, ft, opts);

        tau_fast(i) = (-1/fitresult.b);
        tau_slow(i) = (-1/fitresult.d);
    end

%% Format results tables in workspace

Results = [capacitance, V_hold, I_hold, peak_I_abs, Vhalf, k];
Activation_results = array2table(Results,'VariableNames',{'capacitance (pF)','V_hold (mV)','I_hold (pA)'...
    'peak absolute current (pA)','Vhalf activation (mV)','slope factor k activation'});

Results_per_step = [transpose(V_steps), transpose(peak_I), transpose(peak_I_density),...
    transpose(time_to_peak), transpose(norm_G), transpose(tau_fast), transpose(tau_slow)];
Activation_results_per_step = array2table(Results_per_step,'VariableNames',{'Voltage step (mV)',...
    'PeakI (pA)','PeakI_Density (pA per pF)','TimeToPeak (ms)','Normalized conductance','TauFast (ms)','TauSlow (ms)'});
