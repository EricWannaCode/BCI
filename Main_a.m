close all, clear all, clc

%%

load('EEG.MAT')
fs = eeg.rate; eeg_data = eeg.data; epoch_num = 10;
eeg_data = eeg_data(:,fs*10+1:fs*30); 
FPz_power = zeros(epoch_num,1); Pz_power = zeros(epoch_num,1);
Fz_power = zeros(epoch_num,1); Oz_power = zeros(epoch_num,1); Cz_power = zeros(epoch_num,1);

Pz_relative_power = zeros(epoch_num,5); Cz_relative_power = zeros(epoch_num,1);
Oz_relative_power = zeros(epoch_num,5); Fz_relative_power = zeros(epoch_num,5);
FPz_relative_power = zeros(epoch_num,5);
npoints = 2*fs; fund_frequency = 1/2;
%% select channel
FPz = eeg_data(2,:); Fz = eeg_data(6,:);
Cz = eeg_data(11,:); Pz = eeg_data(16,:);
Oz = eeg_data(20,:); 


FPz = reshape(FPz,npoints,10).';Fz = reshape(Fz,npoints,10).';
Cz = reshape(Cz,npoints,10).';Pz = reshape(Pz,npoints,10).';
Oz = reshape(Oz,10,npoints);
data_all = cat(3,FPz,Fz,Cz,Pz,Oz);
frequency = fs*(0:(npoints/2))./npoints;
names = {'FPz','Fz','Cz','Pz','Oz'};
for ch = 1:5
    data = squeeze(data_all(:,:,ch));
    data_delta = zeros(npoints,10); data_beta = zeros(npoints,10);
    data_theta = zeros(npoints,10); data_alpha = zeros(npoints,10);
    chname = names{ch};
    for e = 1:epoch_num
        data1 =data(e,:).';
        data_fq = fft(data1); data_frequency = data_fq(1:length(data1)/2+1,:);
        data_frequency(2:end-1) = 2*data_frequency(2:end-1);  
        [delta_band, alpha_band, beta_band, theta_band] = select_band(data_frequency,frequency);
        data_delta(:,e) = ifft(delta_band,npoints,'symmetric'); data_alpha(:,e) = ifft(alpha_band,npoints,'symmetric'); 
        data_beta(:,e) = ifft(beta_band,npoints,'symmetric'); data_theta(:,e) = ifft(theta_band,npoints,'symmetric'); 
        
    end
        [ab_power,relative_power]= calculate_bandpower(data_delta,data_theta,data_alpha,data_beta,data.',chname);
    
    
    
end



%% sub-bands
for e = 1:epoch_num
data = [FPz(e,:);Fz(e,:);Cz(e,:);Pz(e,:);Oz(e,:)].';
% data_frequency = abs(fftshift(fft(data,npoints)./length(data))).';


data_fq = fft(data); data_frequency = data_fq(1:length(data)/2+1,:);
data_frequency(2:end-1) = 2*data_frequency(2:end-1);  

frequency = fs*(0:(npoints/2))./npoints;
% frequency = frequency(frequency>=0);
[delta_band, alpha_band, beta_band, theta_band] = select_band(data_frequency,frequency);
data_delta = ifft(delta_band,npoints,'symmetric'); data_alpha = ifft(alpha_band,npoints,'symmetric'); 
data_beta = ifft(beta_band,npoints,'symmetric'); data_theta = ifft(theta_band,npoints,'symmetric'); 


    if e == 1
    figure,
    plot(frequency,abs(data_frequency))
    name = {'FPz','Fz','Cz','Pz','Oz'};
    figure,
    for j = 1:5
    subplot(5,5,j)
    plot((data_delta(:,j)))
    title([ name(j) ' \delta band'])
    if j == 1, ylabel('amplitude'); end
    
    subplot(5,5,5+j)
    plot((data_theta(:,j)))
    title([ name(j) ' \theta band'])
    if j == 1, ylabel('amplitude'); end
    
    subplot(5,5,10+j)
    plot((data_alpha(:,j)))
    title([ name(j) ' \alpha band'])
    if j == 1, ylabel('amplitude'); end
    
    subplot(5,5,15+j)
    plot((data_beta(:,j)))
    title([ name(j) ' \beta band'])
    if j == 1, ylabel('amplitude'); end
        
    subplot(5,5,20+j)
    plot(data(:,j))
    hold on
    plot((data_alpha(:,j))+(data_theta(:,j))+(data_beta(:,j))+(data_delta(:,j)))
    xlabel('sample point')
    if j == 1, ylabel('amplitude'); end
    
    end
    legend('data original','data summation');legend boxoff
%     plot(frequency,data_frequency(:,frequency>=0));

    end
end


%% subfunction 
function [delta_band, alpha_band, beta_band, theta_band] = select_band(data_frequency,frequency_sequence)
    % alpha:8-14Hz, beta:14-30Hz, delta:0-4Hz,theta:4-8Hz
    delta_band = zeros(size(data_frequency)); theta_band = zeros(size(data_frequency)); beta_band = zeros(size(data_frequency)); alpha_band = zeros(size(data_frequency));
    delta_band(0<=frequency_sequence & frequency_sequence<=4,:) = data_frequency(0<=frequency_sequence & frequency_sequence<=4,:);
    alpha_band(8<frequency_sequence & frequency_sequence<=14,:) = data_frequency(8<frequency_sequence & frequency_sequence<=14,:);
    beta_band(14<frequency_sequence & frequency_sequence<=30,:) = data_frequency(14<frequency_sequence & frequency_sequence<=30,:);
    theta_band(4<frequency_sequence & frequency_sequence<=8,:) = data_frequency(4<frequency_sequence & frequency_sequence<=8,:);
end