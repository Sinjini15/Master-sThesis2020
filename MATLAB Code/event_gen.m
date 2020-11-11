%%
clear all;
clc;

%===== Defining variables =====%

fs = 250;
gain = 800.6597;
base = 16;
infant_num = 1;

%===== Get annotations =====%
record_name = strcat('rawdata/infant1_ecg');
ann = rdann(record_name, 'atr');

%%
%===== Generate the events ===%


for i=1:length(ann)
    ind = ann(i);
    N0 = ind-5000;
    N = ind+2500;
    [signal,Fs,tm]=rdsamp(record_name,[],N,N0,4);
    tm = double(tm);
    signal = double(signal);
    ecg_corrected = double((signal - base)/gain);
    med1 = medfilt1(ecg_corrected, 0.2*fs);
    med2 = medfilt1(med1, 0.6*fs);
    ecg_baseline = ecg_corrected - med2;
    ecg_baseline = double(ecg_baseline);
%     [val,j] = max(ecg_baseline);
%     ecg_baseline(j) = 0;
    [qrs_amp_raw,qrs_i_raw,delay,ecg_filtered] = pan_tompkin(ecg_baseline,fs,0);
    data_out = zeros(length(tm),2);
    data_out(:,1) = tm;
    data_out(:,2) = ecg_filtered;
    peak_data = zeros(length(qrs_i_raw),2);
    peak_data(:,1) = qrs_i_raw';
    for k=1:length(qrs_i_raw')
        tm_val = qrs_i_raw(k);
        peak_data(k,2) = tm(tm_val);
    end
    peak_data(:,3) = qrs_amp_raw';
%     filename1 = strcat("/Users/sinjinim/Documents/ASU Semesters/Thesis Sem/My code/event_data/infant",string(infant_num),"/infant",string(infant_num),"_event_",num2str(i),".csv");
%     writematrix(data_out, filename1);
%     filename3 = strcat("/Users/sinjinim/Documents/ASU Semesters/Thesis Sem/My code/processed_data/infant",string(infant_num),"/infant",string(infant_num),"_event_",num2str(i),"_rpeaks.csv");
%     writematrix(peak_data, filename3);
end 

%Generating the .mat files
save event_gen.mat data_out peak_data ann

% filename2 = strcat("/Users/sinjinim/Documents/ASU Semesters/Thesis Sem/My code/event_data/infant", string(infant_num),"/infant",string(infant_num),"_ann.csv");
% writematrix(ann,filename2);




% %--- Plotting the raw ECG signal ----%
% 
% figure;
% plot(signal(1:1000));
% ylabel('Amplitude (mV)');
% xlabel('Indices');
% title('Raw ECG');
% % 
% %--- Plotting baseline wander removed ----%
% 
% figure;
% plot(ecg_baseline(2:1000));
% ylabel('Amplitude (mV)');
% xlabel('Indices');
% title('Processed ECG');
% % 
% % 
% % 
% %--- Pan Tompkin output ----%
% 
% figure
% plot(ecg_filtered(2:2000),'o','MarkerSize', 10,'MarkerIndices',qrs_i_raw(1:17),'MarkerFaceColor','red');
% plot(ecg_filtered(2:2000), 'LineWidth',1.5);
% hold on;
% plot(qrs_i_raw(1:17),qrs_amp_raw(1:17), 'o');
% xlabel("Indices", 'FontSize', 18);
% ylabel("Amplitude (mV)", 'FontSize', 18);
% title('Pan Tompkin output with detected peaks', 'FontSize', 28);
% 
% 
% %----- Plotting to see baseline removal ----%
% 
% figure;
% %subplot (2,1,1);
% plot(signal,'LineWidth',1.5);
% ylabel('Amplitude (mV)','FontSize', 18);
% xlabel('Indices', 'FontSize', 18);
% title('Raw ECG','FontSize', 28);
% %subplot(2,1,2);
% 
% figure;
% plot(ecg_baseline);
% ylabel('Amplitude (mV)', 'FontSize', 18);
% xlabel('Indices', 'FontSize', 18);
% title('ECG with baseline wander removed', 'FontSize', 28);
% 
