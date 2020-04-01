function rtfTemp = generate_template(dpp,freRan,MP)
%% Input
%
% dpp: direct-path propagation, could be i) direct-path impulse response,
% such as HRIR for dummy head case, with the size of (No. of samples) x 
% (No. of microphones) x (No. of candidate locations), and ii) TDOA between
% microphones, with the size of (No. of microphone pairs) x (No. of candidate 
% locations). For both cases, the sampling rate is 16 kHz.
% MP: microphone pairs used for localization, with the size of 2 x (No. of microphone pairs)
% freRan: frequency range used for SSL. It is denoted in frequency bin, within the range of [2 ftLen/2].
%
%% Output
%
% rtfTemp: with size of (No. of frequencies) x (No. of microphone pairs) x (No. of candidate locations)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 16000;
ftLen = 256;                         % Length of Fourier transform

if nargin<2
    freRan = 2:ftLen/4;              % Frequency range up to 4 kHz is used for SSL
end
freNum = length(freRan);

if ndims(dpp) == 3
    [~,micNum,candNum] = size(dpp);
    
    if nargin<3
        MP = [];
        for m1 = 1:micNum-1
            for m2 = m1+1:micNum
                MP = [MP [m1;m2]];
            end
        end
    end
    mpNum = size(MP,2);                      
    
    TF = fft(dpp,ftLen);
    TF = TF(freRan,:,:);
    
    rtfTemp = zeros(freNum,mpNum,candNum);
    for mp = 1:mpNum
        m1 = MP(1,mp);
        m2 = MP(2,mp);
        rtfTemp(:,mp,:) = TF(:,m2,:)./TF(:,m1,:);
    end    
    rtfTemp = rtfTemp./(abs(rtfTemp)+1);  % normalize level to [0,1]
else
    [mpNum,candNum] = size(dpp);    
    rtfTemp = 0.5*exp(-1i*bsxfun(@times,repmat(2*pi*fs*freRan'/ftLen,[1,mpNum,candNum]), reshape(dpp,[1,mpNum,candNum])));  
end

