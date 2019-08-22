function OUT = processImpulseResponse(imp_param)
% OUT is a structure that contains the following variables:
% OUT.DeltaL_meas, the level difference between the two microphones as
% measured in the experiment. This is the key data to be used in MSMNT
% method. This is a function of frequency band "OUT.foct".
%
% OUT.foct, the frequency band considered for the computation of the level
% difference, i.e. "OUT.DeltaL_meas".
%
% OUT.config_dim, the variable containing the physical dimensions of the
% experimental setup. It reads from the "config.csv" that is available in
% the same folder of the data. "config.csv" contains 4 rows and multiple
% columns, depending on the number of configurations. The rows corresponds
% to:
% Row 1: height of source, meter
% Row 2: height of microphone 1, meter
% Row 3: height of microphone 2, meter
% Row 4: distance between source and microphone
%
%
% INPUT: imp_param contains in the following parameters:
% imp_param.datafolder, where are all the filed located. The files should
% have the following syntax, e.g. C1M11.wav, indicating 
% C[config_number]M[microphone_number][measurement_index], where C and M
% stands for Configuration and Microphone. For the 4th measurement of the
% 2nd configuration on 2nd microphone, the file willl be C2M24.wav.
%
% imp_param.n_layer, number of layers to be computed. For just a
% semi-infinite layer, e.g. concrete, n_layer = 0. For a surface that is 
% on top of a hard-backed layer, n_layer = 1. For a dual layers, e.g. green
% roof, n_layer = 2. See Chang Liu's thesis Pg 10 Figure 1.5 for an
% illustrative example.
% thesis link: https://research.tue.nl/nl/publications/in-situ-characterization-of-the-acoustic-impedance-of-vegetated-r
%
% imp_param.guess_thickness, 
% an estimation of the total thickness of all the layers on the surface. 
% This does not have to be very accurate as it is only used to truncate 
% the impulse responses.
%
% imp_param.n_meas, 
% number of measurements taken. This corresponds to the maximum number 
% of the last digit in the date files.
%
% imp_param.fs, 
% sampling frequency during data acquisition phase.
%
% imp_param.n_config, 
% number of unique configurations performed. 
% Recommendaded minimum 3 for n_layer = 1; minimum 6 for n_layer = 2.
%
% imp_param.flow, 
% lowest considered frequency. Lowest as recommended by the Nordtest 
% method is 200 Hz. The script considers only the frequency bands that 
% are bounded by both "flow" and "fhigh" (inclusive).
%
% imp_param.fhigh,
% highest considered frequency. Highest as recommended by the Nordtest 
% method is 2,500 Hz. The script considers only the frequency bands 
% that are bounded by both "flow" and "fhigh" (inclusive).
%
% V1: 22/8/19, J.J. Tan.

datafolder = imp_param.datafolder;
n_meas     = imp_param.n_meas;
fs         = imp_param.fs;
n_config   = imp_param.n_config;
flow       = imp_param.flow;
fhigh      = imp_param.fhigh;
n_layer    = imp_param.n_layer;

if n_layer ~= 0
    guess_thickness = imp_param.guess_thickness;
else
    guess_thickness = 0;
end

% we calculate an estimated path of travel of the sound waves
config_dim = load(fullfile(datafolder,'config.csv'));
hs  = config_dim(1,:) - guess_thickness;
hr1 = config_dim(2,:) - guess_thickness;
hr2 = config_dim(3,:) - guess_thickness;
dsr = config_dim(4,:);

R1{1}      = sqrt( ( hs -  hr1).^2 +  dsr.^2 ); %direct path of mic 1
R2{1}      = sqrt( ( hs +  hr1).^2 +  dsr.^2 ); %reflected path of mic 1
R1{2}      = sqrt( ( hs -  hr2).^2 +  dsr.^2 ); %direct path of mic 2
R2{2}      = sqrt( ( hs +  hr2).^2 +  dsr.^2 ); %reflected path of mic 1
diffdistM1 = R2{1} - R1{1};
diffdistM2 = R2{2} - R1{2};

for kk = 1:n_config 

    %here we extract the impulse response recordings from .wav files
    for i = 1:n_meas
        filename_mic1 = ['C',num2str(kk),'M1',num2str(i),'.wav'];
        filename_mic2 = ['C',num2str(kk),'M2',num2str(i),'.wav'];

        [M1(:,i),Fs] = audioread(fullfile(datafolder,filename_mic1));
        [M2(:,i),Fs] = audioread(fullfile(datafolder,filename_mic2));

    end

    %then we have to trim the signal to include only first incident and first
    %reflection pulses only
    %by choice due to ease, we will employ a linear fade out function...
    %linear fade out looks somewhat like...  1 _______
    %                                                 \__
    %                                        0
    %1 at the part we are interested in, then a linear slope that gradually 
    % goes to 0.
    
    %NOTE: current implementation is a bit clunky but for what it is worth,
    %it works. this will be revisited in the future.
    
    % take the maximum value of the signal - this is the initial sound wave 
    % that arrives at the microphone directly
    [~,ind_firstpulseM1] = max(M1);
    [~,ind_firstpulseM2] = max(M2);
    

    % we then calculate the number of points between the direct path and
    % reflected path based on the extra travel distance of reflected path.
    % ind_firstreflect = extra travel distance / speed of sound * sampling freq;
    ind_firstreflectM1 = ceil(diffdistM1(kk)/343*fs); %relative to first pulse
    ind_firstreflectM2 = ceil(diffdistM2(kk)/343*fs); %relative to first pulse
    
    % the cut point is:
    % where the first pulse is + distance to second pulse + 0.5ms (fs/2e3)
    % the 0.5ms signal more is to ensure the information of second pulse is
    % properly included.
    i_cutpointM1 = ind_firstpulseM1 + ind_firstreflectM1 + fs/2e3; 
    i_cutpointM2 = ind_firstpulseM2 + ind_firstreflectM2 + fs/2e3; 

    % by choice, we define the linear slope to be 0.5ms long, i.e. fs/2e3
    % thus, complete silence will start after the cutpoint + fade_length,
    % as defined by i_silenceM1 and i_silenceM2.
    fade_length = fs/2e3;
    i_silenceM1 = i_cutpointM1 + fade_length;
    i_silenceM2 = i_cutpointM2 + fade_length;
    
    % however, not all i_silenceM1/M2 are the same (compared to the
    % different measurement set, i.e. [measurement_index] in the .wav file)
    % for this, we define the number of zero padding needed ...
    n_zeropadM1 = max(i_cutpointM1) - i_cutpointM1;
    n_zeropadM2 = max(i_cutpointM2) - i_cutpointM2;
    
    % and create a new matrix of the correct length
    M1t = zeros(max(i_cutpointM1)+1 + fade_length,n_meas);
    M2t = zeros(max(i_cutpointM2)+1 + fade_length,n_meas);

    % and run a loop through all mesaurement_index
    for ww = 1:1:n_meas
        M1t(:,ww) = [M1(1:i_cutpointM1(ww),ww); ... %all the important data points
            M1(i_cutpointM1(ww)+1:i_silenceM1(ww),ww).*linspace(1,0,fade_length)'; ... %all the fading data points
            zeros(n_zeropadM1(ww)+1,1)]; %all the padded 0 data points, if necessary
        M2t(:,ww) = [M2(1:i_cutpointM2(ww),ww); ... %all the important data points
            M2(i_cutpointM2(ww)+1:i_silenceM2(ww),ww).*linspace(1,0,fade_length)'; ... %all the fading data points
            zeros(n_zeropadM2(ww)+1,1)]; %all the padded 0 data points, if necessary
    end
    
    
    % then, we can determine the level difference.
    % first, we perform FFT on our trimmed time signals
    M1f = fft(M1t,fs); %second argument pad output to length fs
    M2f = fft(M2t,fs); 

    % then we take the mean of the spectra
    L_meas1 = mean(abs(M1f),2);
    L_meas2 = mean(abs(M2f),2);

    % and we compute the octave band values
    [foct,Lp_meas1] = compute_octave_band_SPL(1:1:fs,L_meas1,flow,fhigh,'1/3');
    [~,Lp_meas2] = compute_octave_band_SPL(1:1:fs,L_meas2,flow,fhigh,'1/3');

    % and we take the difference, because we are looking for 
    % level _difference_
    DeltaL_meas(:,kk) =  Lp_meas1 - Lp_meas2;

end

OUT.DeltaL_meas = DeltaL_meas;
OUT.foct        = foct;
OUT.config_dim  = config_dim;

