clear all;
close all;
clc;

% Load the ECG data
ECG = load('samin_ECG-L05');
ECG_1 = ECG.data(:,1);  % Channel 1 (ECG signal)

% Assume the sampling rate is given
Fs = 1011.35;  % Sampling rate in Hz

% Define the time vector based on the number of samples and sampling rate
num_samples = length(ECG_1);
ts = (0:num_samples-1) / Fs;  % Time vector in seconds

% Stage timings in seconds
stages = [0, 19.20;   % Supine
          19.21, 38.40;   % Seated
          38.41, 70.40;   % Deep Breathing
          70.41, 140.80]; % After Exercise

stage_names = {'Supine', 'Seated', 'Deep Breathing', 'After Exercise'};

% Bandpass filter for ECG signal (0.5 Hz to 40 Hz)
low_cutoff = 0.5;   % Lower bound of heart rate (0.5 Hz)
high_cutoff = 40;   % Upper bound of heart rate (40 Hz)
[b, a] = butter(2, [low_cutoff, high_cutoff] / (Fs / 2), 'bandpass');  % 2nd-order Butterworth filter

% Loop through each stage
for i = 1:size(stages, 1)
    
    % Convert stage times to indices
    start_idx = max(1, round(stages(i, 1) * Fs));  % Ensure start index is at least 1
    end_idx = min(num_samples, round(stages(i, 2) * Fs));  % Ensure end index doesn't exceed array size
    
    ECG_stage = ECG_1(start_idx:end_idx);
    ts_stage = ts(start_idx:end_idx);

    % Filter the ECG signal
    filtered_ECG_stage = filtfilt(b, a, ECG_stage);

    % Plot original ECG data for the current stage
    figure;
    subplot(2, 1, 1);
    plot(ts_stage, ECG_stage);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(['Original ECG Signal - ', stage_names{i}]);
    grid on;
    
    % Plot filtered ECG data for the current stage
    subplot(2, 1, 2);
    plot(ts_stage, filtered_ECG_stage);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(['Filtered ECG Signal - ', stage_names{i}]);
    grid on;
    
    % Adjust R-peak detection parameters based on the stage
    if strcmp(stage_names{i}, 'After Exercise')
        [~, r_locs] = findpeaks(filtered_ECG_stage, 'MinPeakHeight', 0.04, 'MinPeakDistance', 0.4 * Fs);  % More sensitive for high HR
    elseif strcmp(stage_names{i}, 'Deep Breathing')
        [~, r_locs] = findpeaks(filtered_ECG_stage, 'MinPeakHeight', 0.004, 'MinPeakDistance', 0.44 * Fs);
    else
        [~, r_locs] = findpeaks(filtered_ECG_stage, 'MinPeakHeight', 0.05, 'MinPeakDistance', 0.37 * Fs);  % Standard peak detection
    end

    % Mark R-peaks on filtered signal
    figure;
    plot(ts_stage, filtered_ECG_stage);
    hold on;
    plot(ts_stage(r_locs), filtered_ECG_stage(r_locs), 'ro');  % Mark R-peaks with red circles
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(['Filtered ECG Signal with Detected R-peaks - ', stage_names{i}]);
    grid on;

    % Calculate BPM for the current stage
    if length(r_locs) > 1
        rr_intervals = diff(r_locs) / Fs;  % Time intervals between R-peaks in seconds
        average_rr = mean(rr_intervals);   % Average time between R-peaks
        bpm = 60 / average_rr;             % Convert to beats per minute
    else
        bpm = NaN;  % Not enough peaks detected
    end

    % Display BPM result
    disp([stage_names{i}, ' - Estimated Heart Rate: ', num2str(bpm), ' BPM']);
end
