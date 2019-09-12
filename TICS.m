% Created by Simon Sehayek
% This function is meant to take an image series and calculate its
% TICS autocorrelation function.
%
% INPUT VARIABLES
%
%   imgser: image series
%
% OUTPUT VARIABLES
%
%   r_tau_avg: averaged TICS autocorrelation function of imgser
%
% VARARGIN
%
% 'bkgdCxn','bkgd'
% -------------------------------------------------------------------------
% 0 (default) | array
% -------------------------------------------------------------------------
% This is only relevant if using 'ticsNorm' option set to 1. The background
% subtracted image series is recommended for computing the mean used in the
% normalization, in this case. Supply background subtracted image series if
% using this option. By default, the mean of imgser is used in the
% normalization.
%
% 'computationMethod','compMethod','method','useMethod'
% -------------------------------------------------------------------------
% 'WKT','WKTMethod' (default) | 'explicit'
% -------------------------------------------------------------------------
% Method for computing TICS autocorrelation. By default, Wiener-Khinchin is
% used (recommended). The other option explicitly computates the
% autocorrelation through its definition.
%
% 'maxLag'
% -------------------------------------------------------------------------
% n_frames-1 (default) | positive integer (value does not exceed
% n_frames-1)
% -------------------------------------------------------------------------
% Maximum number of time lags to compute in autocorrelation. Only relevant
% when 'computationMethod' is set to 'explicit'. By default, all time lags
% are computed.
%
% 'meanType'
% -------------------------------------------------------------------------
% 'spatial' (default) | 'temporal'
% -------------------------------------------------------------------------
% Type of intensity fluctuations to autocorrelate. By default, the spatial
% mean is subtracted from intensity values to compute fluctuations.
%
% 'parallel','useParallel'
% -------------------------------------------------------------------------
% 0 (default) | 1
% -------------------------------------------------------------------------
% Boolean for parallelization of TICS autocorrelation computation. This
% option should be used for image series with many pixels, especially if
% MATLAB runs out of memory.
%
% 'ticsNorm'
% -------------------------------------------------------------------------
% 0 (default) | 1
% -------------------------------------------------------------------------
% Boolean for using traditional TICS normalization.[1] By default, this is
% not done.
%
% 'windowAvg','winAvg','movAvg','windowMean','winMean','movMean','window'
% -------------------------------------------------------------------------
% 0 (default) | positive integer
% -------------------------------------------------------------------------
% Parameter k or [kb,kf] which goes into movmean. By default, no windowed
% mean is performed. We reccomend choosing k which is slightly bigger than
% the PSF size. This is helpful in reducing spatial non-uniformity in the
% background intensity.
%
% REFERENCES
%
% [1] Kolin, D. L.; Wiseman, P. W. Cell Biochemistry and Biophysics 2007,
%     49, 141-164.
%
function [r_tau_avg] = TICS(imgser,varargin)

% default values for varargin options
tics_norm = 0;
mean_type = 'spatial';
win_mean = 0;
comp_method = 'WKT';
use_parallel = 0;
%
bkgd_cxn = 0;
J_sub = [];
%

size_x = size(imgser,2);
size_y = size(imgser,1);
n_frames = size(imgser,3);
max_lag = n_frames - 1;

for i = 1:2:length(varargin)
    if any(strcmpi(varargin{i},{'ticsNorm'}))
        if isnumeric(varargin{i+1}) && any(varargin{i+1} == [0,1])
            tics_norm = varargin{i+1};
        else
            warning(['Unknown option for ''',varargin{i},...
                ''', using default options.'])
        end
    elseif strcmpi(varargin{i},'meanType')
        if strcmpi(varargin{i+1},'spatial')
            mean_type = 'spatial';
        elseif strcmpi(varargin{i+1},'temporal')
            mean_type = 'temporal';
        else
            error(['Invalid mean type specified; specify either',...
                ' ''spatial'' or ''temporal''.'])
        end
    elseif any(strcmpi(varargin{i},{'windowAvg','winAvg','movAvg',...
            'windowMean','winMean','movMean','window'}))
        % option to compute windowed mean for subtraction and normalization
        % (if applicable). Enter a positive integer for an equal back/forth
        % windowing, enter length 2 vector for unequal back/forth
        % windowing, [kb kf] (see movmean doc). In the case of spatial
        % mean type, the windowing will be the same in both dimensions (can
        % be generalized if needed)
        % another generalization can include windowed spatio-temporal means
        if isnumeric(varargin{i+1}) && varargin{i+1} > 0 && ...
                length(varargin{i+1})<=2
            % set windowing boolean to 1
            win_mean = 1;
            % window size specification
            win_k = varargin{i+1};
        else
            error(['Invalid mean type specified; specify either',...
                ' ''spatial'' or ''temporal''.'])
        end
    elseif any(strcmpi(varargin{i},{'computationMethod','compMethod',...
            'method','useMethod'}))
        if any(strcmpi(varargin{i+1},{'explicit'}))
            comp_method = 'explicit';
        elseif any(strcmpi(varargin{i+1},{'WKT','WKTMethod'}))
            comp_method = 'WKT';
        else
            warning(['Unknown option for ''',varargin{i},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{i},{'parallel','useParallel'}))
        % use this option if MATLAB runs out of memory
        if any(varargin{i+1}==[0,1])
            use_parallel = varargin{i+1};
        else
            warning(['Unknown option for ''',varargin{i},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{i},{'bkgdCxn','bkgd'}))
        if isnumeric(varargin{i+1}) && isequal(size(varargin{i+1},3),n_frames)
            bkgd_cxn = 1;
            J_sub = varargin{i+1};
        else
            warning(['Unknown option for ''',varargin{i},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{i},{'maxLag'}))
        if isnumeric(varargin{i+1}) && varargin{i+1} >= 0
            max_lag = varargin{i+1};
            disp(['max lag set to ',num2str(max_lag),' out of ',num2str(n_frames-1)])
        else
            warning(['Unknown option for ''',varargin{i},...
                ''', using default options.'])
        end
    else
        warning(['unknown varargin input ''',varargin{i},'''.'])
    end
end

% compute mean as specified by mean_type and win_mean
if win_mean
    % windowed mean
    if strcmp(mean_type,'spatial')
        mean_im = movmean(movmean(imgser,win_k,1,'omitnan'),win_k,2,'omitnan');
        if bkgd_cxn
            mean_im_sub = movmean(movmean(J_sub,win_k,1,'omitnan'),win_k,2,'omitnan');
        end
    else
        mean_im = movmean(imgser,win_k,3);
        if bkgd_cxn
            mean_im_sub = movmean(J_sub,k,3);
        end
    end
else
    % standard mean
    if strcmp(mean_type,'spatial')
        mean_im = repmat(mean(mean(imgser,1,'omitnan'),2,'omitnan'),...
            [size_y,size_x,1]);
        if bkgd_cxn
            mean_im_sub = repmat(mean(mean(J_sub,'omitnan'),'omitnan'),...
                [size_y,size_x]);
        end
    else
        mean_im = repmat(mean(imgser,3),[1,1,n_frames]);
        if bkgd_cxn
            mean_im_sub = repmat(mean(J_sub,3),[1,1,n_frames]);
        end
    end
end

% compute image intensity fluctuations according to tics_norm and bkgd_cxn
if tics_norm && bkgd_cxn
    fluct_imgser = (imgser-mean_im)./mean_im_sub;
elseif tics_norm && ~bkgd_cxn
    fluct_imgser = (imgser-mean_im)./mean_im;
else
    fluct_imgser = imgser-mean_im;
end

% clear imgser as it is no longer needed
clear imgser
if strcmp(comp_method,'WKT')
    % find next power of 2 before FFT
    if floor(log2(n_frames)) == log2(n_frames)
        T_pad = 2*n_frames;
    else
        pow2 = nextpow2(2*n_frames);
        T_pad = 2^pow2;
    end
    %
    
    % FFT of intensity fluctuations
    if ~use_parallel
        fourierFluct = fft(fluct_imgser,T_pad,3);
        r_tau = ifft(fourierFluct.*conj(fourierFluct),[],3);
    else
        parpool
        
        r_tau = zeros(size_y,1,T_pad);
        parfor y = 1:size_y
            % fourier transform all pixels in row y
            fourierFluct_y = fft(fluct_imgser(y,:,:),T_pad,3);
            % inverse fourier transform all pixels in row y
            r_tau_y = ifft(fourierFluct_y.*conj(fourierFluct_y),[],3);
            % take mean along all pixels in row y
            r_tau(y,:,:) = mean(r_tau_y,2);
        end
        delete(gcp)
        
        r_tau = mean(r_tau,1);
    end
    
    r_tau(:,:,n_frames+1:end) = [];
    for tau = 0:n_frames-1
        r_tau(:,:,tau+1) = 1./(n_frames-tau).*r_tau(:,:,tau+1);
    end
    
else
    r_tau = zeros(size_y,size_x,max_lag+1);
    tic
    for tau = 0:max_lag
        r_tau(:,:,tau+1) = 1/(n_frames-tau)*dot(fluct_imgser(:,:,1:n_frames-tau),fluct_imgser(:,:,tau+1:n_frames),3);
        progBar(tau,max_lag,'time',toc);
    end
end

% taking absolute value biases the autocorrelation when it is close to 0
%
% r_tau_avg = abs(mean(mean(r_tau,1),2));
r_tau_avg = mean(mean(r_tau,1,'omitnan'),2,'omitnan');
r_tau_avg = squeeze(r_tau_avg);