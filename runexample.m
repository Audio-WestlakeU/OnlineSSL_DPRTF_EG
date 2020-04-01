clear 

% The sampling rate of signals/RIRs is set to 16 kHz, frame length for STFT is set to 256 samples (16 ms)
fs = 16000;
ftLen = 256;               

MP = [1 1 1 2 2 3
      2 3 4 3 4 4];         % all the six microphone pairs are used

%% Kinovis-MST example

% load HRIR
load('data/HRIR_NAO_48kHz_-175:5:180_degrees.mat');

freRan = [26:31 34:107];    % frequency bins used for Kinovis-MST example, which are not too noisy

% generate template  
rtfTemp = generate_template(HRIR_NAO,freRan,MP);

% read audio 
[y,rfs] = audioread('data/Kinovis-MST-example.wav');
micNum = size(y,2);
x = [];
for mic = 1:micNum
    x(:,mic) = resample(y(:,mic),fs,rfs);
end

% localization
[GMMWeight,Peaks] = OnlineSSL_DPRTF_EG(x,rtfTemp,freRan,MP);

% resultplot

load('data/Kinovis-MST-speakerPosition.mat') 
load('data/Kinovis-MST-speakerVAD.mat')
speakerPositionVAD = speakerPosition.*(speakerVAD);
speakerPositionVAD = speakerPositionVAD + (speakerPositionVAD==0)*200;

fraNum = size(speakerPosition,1);

figure;

subplot(311);hold;
plot(speakerPositionVAD,'.')
axis([1 fraNum -175 180])
set(gca,'xtick',1250:1250:fraNum,'xticklabels',10:10:fraNum/125,'ytick',-120:60:120,'yticklabels',-120:60:120,'FontSize',12,'box','on')
title('ground truth')

subplot(312);imagesc(GMMWeight(end:-1:1,:))
set(gca,'xtick',1250:1250:fraNum,'xticklabels',10:10:fraNum/125,'ytick',12:12:60,'yticklabels',-120:60:120,'FontSize',12,'box','on')
title('GMM Weight')
ylabel('Azimuth (degrees)')

subplot(313);imagesc(Peaks(end:-1:1,:))
set(gca,'xtick',1250:1250:fraNum,'xticklabels',10:10:fraNum/125,'ytick',12:12:60,'yticklabels',-120:60:120,'FontSize',12,'box','on')
title('Peak selection')
xlabel('Time (s)')


%% LOCATA example

% we use four (out of twelve) microphones
micPosition = [-0.031  0.023 0.042;
               -0.036 -0.027 0.038;
                0.034 -0.030 0.037;
                0.035  0.025 0.039];

% candidate locations
AZI = (-pi+pi/36:pi/36:pi)';      % azimuth -175:5:180
ELE = pi*ones(size(AZI))/2;        
Ran = 10;
canPosition = Ran*[-sin(ELE).*sin(AZI),sin(ELE).*cos(AZI),cos(ELE)]+0.04;

% TDOA
TDOA = compute_TDOA(micPosition,canPosition,MP);

% generate template
freRan = 2:65;  
rtfTemp = generate_template(TDOA,freRan);

% read audio, four channels ([5 8 11 12]) of LOCATA Robot recording
[y,rfs] = audioread('data/LOCATA-dev-task6-rec3.wav');
micNum = size(y,2);
x = [];
for mic = 1:micNum
    x(:,mic) = resample(y(:,mic),fs,rfs);
end


% localization
[GMMWeight,Peaks] = OnlineSSL_DPRTF_EG(x,rtfTemp,freRan,MP);

% result plot

load('data/LOCATA-speakerPosition.mat')
load('data/LOCATA-speakerVAD.mat')
speakerPositionVAD = speakerPosition.*(speakerVAD);
speakerPositionVAD = speakerPositionVAD + (speakerPositionVAD==0)*200;

fraNum_gt = size(speakerPosition,1);
fraNum = size(GMMWeight,2);

figure;

subplot(311);hold;
plot(speakerPositionVAD,'.')
axis([1 fraNum_gt -175 180])
set(gca,'xtick',1200:1200:fraNum_gt,'xticklabels',10:10:fraNum_gt/120,'ytick',-120:60:120,'yticklabels',-120:60:120,'FontSize',12,'box','on')
title('ground truth')

subplot(312);imagesc(GMMWeight(end:-1:1,:))
set(gca,'xtick',1250:1250:fraNum,'xticklabels',10:10:fraNum/125,'ytick',12:12:60,'yticklabels',-120:60:120,'FontSize',12,'box','on')
title('GMM Weight')
ylabel('Azimuth (degrees)')

subplot(313);imagesc(Peaks(end:-1:1,:))
set(gca,'xtick',1250:1250:fraNum,'xticklabels',10:10:fraNum/125,'ytick',12:12:60,'yticklabels',-120:60:120,'FontSize',12,'box','on')
title('Peak selection')
xlabel('Time (s)')

