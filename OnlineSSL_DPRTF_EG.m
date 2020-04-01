function [GMMWeight,Peaks] = OnlineSSL_DPRTF_EG(soundData,rtfTemp,freRan,MP)
%% Input
% Sound Data : audio signal, with the size of (No. of samples) x (No. of mics), and with the sampling rate of 16 kHz
% rtfTemp: RTF template for candidate locations, with the size of (No. of freqencies) x (No. of microphone pairs) x (No. of candidate locations)
% freRan: frequency range used for SSL. It is denoted in frequency bin, within the range of [2 winLen/2].
% MP: microphone pairs used for localization, with the size of 2 x (No. of microphone pairs)
%
%% Output
% GMMWeight: (No. of candidate locations) x (No. of frames).
% Peaks: (No. of candidate locations) x (No. of frames), with 1 for detected active speakers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an online implementation for RLS DP-RTF estimation and EG multispeaker localization
%
% The methods are described in the paper:
%
% - X. Li, Y. Ban, L. Girin, X. Alameda-Pineda and R. Horaud. Online Localization
% and Tracking of Multiple Moving Speakers in Reverberant Environments. IEEE JSTSP 2019.
%
% and related papers:
%
% - X. Li, L. Girin, R. Horaud and S. Gannot. Estimation of the Direct-Path
% Relative Transfer Function for Supervised Sound-Source Localization. TASLP 2016.
% - X. Li, L. Girin, R. Horaud and S. Gannot. Multiple-Speaker Localization Based on
% Direct-Path Features and Likelihood Maximization with Spatial Sparsity
% Regularization. TASLP 2017.
% - X. Li, B. Mourgue, L. Girin, S. Gannot and R. Horaud. Online Localization
% of Multiple Moving Speakers in Reverberant Environments. SAM 2018.
%
% Author: Xiaofei Li, INRIA Grenoble Rhone-Alpes
% Copyright: Perception Team, INRIA Grenoble Rhone-Alpes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[d1,d2] = size(soundData);
if d1<d2
    soundData = soundData';
end
micNum = size(soundData,2);

%% Microphone pairs
if nargin<4
    MP = [];
    for m1 = 1:micNum-1
        for m2 = m1+1:micNum
            MP = [MP [m1;m2]];
        end
    end
end
mpNum = size(MP,2);

%% STFT parameters
fs = 16000;
winLen = 256;
fraInc = winLen/2;
win=hamming(winLen);

if nargin<3
    freRan = 2:winLen/4;                                 % Frequency range up to 4 kHz
end
freNum = length(freRan);

fraNum = floor((length(soundData)-winLen)/fraInc);       % Number of frames

%% Spectral smoothing and Speech/non-speech classification parameters
smoLen = ceil(200*16/fraInc);                     % Spectral smoothing history length, 200 ms
beta = (smoLen-1)/(smoLen+1);                     % Smoothing factor
minSran = 3;                                      % Window length for minimum searching, in senconds
mcmtFn = floor(minSran*fs/fraInc)-smoLen+1;
[sTh,nTh] = MCMT(smoLen,1,mcmtFn);                % Thresholds for speech and noise frames detection

sfreNumTh = 0.05*freNum;                          % If speech frequency is less than this number, skip

%% RLS DP-RTF estimation parameters
% IMPORTANT! CTF length, this setting is suitable for reverberation time around 0.55 s, it approximately equals T60/8. 
% It has a large influence to localization performance, users could adjust this accordingly to reverberation time.
fiLen = 8;                                        

rctfLen = fiLen*micNum;                           % Length of relative CTF vector
fraNumLs = ceil(fiLen/sqrt(micNum-1));            % Number of frames used for relative CTF estimation
lambda = (fraNumLs-1)/(fraNumLs+1);               % forgetting factor

ctTh = 0.5;                                       % Threshold for consistency test
%% EG parameters
candNum = size(rtfTemp,3);                        % Number of candidate locations
noInfVec = ones(candNum,1)/candNum;               % Uniform vector

eta = 0.07;                                       % EG update rate
eta1 = 0.065;

gamma = 0.1;                                      % Entropy regularization weight
spaSmooth = 0.02;                                 % Spatial smoothing factor

%% Peak detection parameters
psTh = 0.04;                                      % IMPORTANT! peak detection threshold, controls the tradeoff between false alarm and miss detection
maxPeak = 5;                                      % maximum number of active speakers

%% Variables in process

ftHis = zeros(freNum,micNum,fiLen);               % STFT coefficients history
sft = zeros(freNum,micNum,fiLen);                 % Smoothed auto-/cross-spectra
nearNoiseFra = zeros(freNum,micNum,fiLen);        % The nearest noise frame
ftPow = zeros(freNum,1);                          % Smoothed STFT power

covInv = zeros(freNum,rctfLen,rctfLen);           % Signal covariance matrix inverse
for fre = 1:freNum
    covInv(fre,:,:) = 1e3*eye(rctfLen);
end

gmmWeight = noInfVec;                             % GMM weight vector in process
GMMWeight = zeros(candNum,fraNum);                % GMM weight for all frames

Peaks = zeros(candNum,fraNum);                    % peaks

%% DP-RTF estimation and EG multispeaker localization, frame by frame

WH = waitbar(0,'Please wait ...');
for t = 1:fraNum
    
    %% signal processing
    
    % Short-time time-domain signal
    xt = soundData((t-1)*fraInc+1:(t-1)*fraInc+winLen,:);
    
    % FFT
    xft = fft(bsxfun(@times,xt,win));
    xft = xft(freRan,:);
    
    % FT coefficients history
    ftHis(:,:,2:fiLen) = ftHis(:,:,1:fiLen-1);
    ftHis(:,:,1) = xft;
    
    % Smoothing
    sft = beta*sft+(1-beta)*bsxfun(@times,ftHis,conj(xft(:,1)));
    
    % The psd of the current frame
    ftPow = beta*ftPow+(1-beta)*sum(abs(xft).^2,2);
    
    % Initialization 
    if t<smoLen
        gmmWeight = (1-eta1)*gmmWeight + eta1*noInfVec;
        GMMWeight(:,t) = gmmWeight;
        continue;        
    elseif t==smoLen
        minPowHis = repmat(ftPow,[1,mcmtFn]);
        gmmWeight = (1-eta1)*gmmWeight + eta1*noInfVec;
        GMMWeight(:,t) = gmmWeight;
        continue;
    end
    
    % Minimum power searching
    minPowHis(:,1:mcmtFn-1) = minPowHis(:,2:mcmtFn);
    minPowHis(:,mcmtFn) = ftPow;
    minPow = min(minPowHis,[],2);   
    
    nearNoiseFra(ftPow<=nTh*minPow,:,:) = sft(ftPow<=nTh*minPow,:,:);    % The nearest noise frame with respect to the current frame
    sftSS = sft-nearNoiseFra;         % Spectral subtraction
    
    % Speech frequency index, skip to next frame if empty
    sInd = ftPow>sTh*minPow;
    sfreNum = sum(sInd);
    if sfreNum < sfreNumTh
        gmmWeight = (1-eta1)*gmmWeight + eta1*noInfVec;
        GMMWeight(:,t) = gmmWeight;
        continue;
    end        
    sftSS = sftSS(sInd,:,:);          % speech frequency bins
    
    %% DP-RTF estimation 
    
    covInvt = covInv(sInd,:,:);     % Signal covariance matrix inversion to be updated
    
    % Signal covariance matrix inversion update based on Sherman-Morrison formula
    for m1 = 1:micNum-1
        for m2 = m1+1:micNum
            xvec = zeros(sfreNum,rctfLen);
            xvec(:,(m1-1)*fiLen+1:m1*fiLen) = squeeze(sftSS(:,m2,:));
            xvec(:,(m2-1)*fiLen+1:m2*fiLen) = -squeeze(sftSS(:,m1,:));
            
            gain = sum(bsxfun(@times,covInvt,reshape(conj(xvec),[sfreNum,1,rctfLen])),3);
            gain = gain./repmat(1+sum(gain.*xvec,2),[1,rctfLen]);
            covInvtDelta = bsxfun(@times,repmat(gain,[1,1,rctfLen]),sum(bsxfun(@times,covInvt,xvec),2));
            covInvt = covInvt-covInvtDelta;
        end
    end
    covInv(sInd,:,:) = covInvt/lambda;
    
    % Extract DP-RTF from the Relative CTF estimates, one column for one
    % reference channel setting, one row for one channel
    dprtfRef = covInvt(:,1:fiLen:end,1:fiLen:end);  
    
    % Consistency test
    dprtf = ones(sfreNum,mpNum);    
    for mp = 1:mpNum
        m1 = MP(1,mp);
        m2 = MP(2,mp);
        a = dprtfRef(:,m2,m1)./dprtfRef(:,m1,m1);
        b = dprtfRef(:,m2,m2)./dprtfRef(:,m1,m2);
        simi = abs(1+a.*conj(b))./sqrt((1+abs(a).^2).*(1+abs(b).^2));
        dprtf(:,mp) = (simi>ctTh).*(a+b)/2;
    end
    dprtf = dprtf./(1+abs(dprtf));   % dprtf estimate, with 0 value for invalid elements   
    
    if sum(sum(dprtf~=0))==0
        gmmWeight = (1-eta1)*gmmWeight + eta1*noInfVec;
        GMMWeight(:,t) = gmmWeight;
        continue;
    end
    
    %% EG multispeaker localization      
    
    loglik = loglik_complexgaussian(dprtf,rtfTemp(sInd,:,:),0.1);  % log-likelihood
    lik = exp(loglik);
    
    dprtf = reshape(dprtf,[sfreNum*mpNum,1]);
    
    lik = reshape(lik,[sfreNum*mpNum,candNum]);
    lik = lik(dprtf~=0,:);         % likelihood for valid dprtf features    
    
    gmmLik = bsxfun(@times,lik,gmmWeight');
    DerLik = -mean(bsxfun(@times,lik,1./(sum(gmmLik,2))),1)';     % derivative of log-likelihood
    
    DerEntropy = -(1+log(gmmWeight));        % derivative of entropy
    
    Der = DerLik+gamma*DerEntropy;           % overall derivative
    
    expDer = exp(-eta*Der);
    expDer(expDer>1e2) = 1e2;                % exponentiated gradient, limited to be less than 100
    
    gmmWeight = gmmWeight.*expDer;           % GMM weight update
    gmmWeight = gmmWeight/sum(gmmWeight);
    
    %  circular/spatial smoothing
    gmmWeightM1 = [gmmWeight(2:end);gmmWeight(1)];
    gmmWeightP1 = [gmmWeight(end);gmmWeight(1:end-1)];
    gmmWeight = (gmmWeight+spaSmooth*(gmmWeightM1+gmmWeightP1))/(1+2*spaSmooth);    
    GMMWeight(:,t) = gmmWeight;
    
    %% Circular peak detection/selection
    peaks = gmmWeight>[gmmWeight(end); gmmWeight(1:end-1)] & gmmWeight>[gmmWeight(2:end); gmmWeight(1)];
    if sum(peaks)>maxPeak
        sortWeight = sort(gmmWeight(peaks),'descend');
        Peaks(:,t) = peaks & gmmWeight>max(psTh,sortWeight(maxPeak+1));
    else
        Peaks(:,t) = peaks & gmmWeight>psTh;
    end
    
    waitbar(t/fraNum)
end
close(WH)