% Written by Simon Sehayek
%
% Generate image series of immobile fluorescent probes undergoing
% blinking and bleaching
%
% INPUT VARIABLES
%
%   sz: number of pixels in each dimension (only square image series are
% supported)
%   T: number of frames
%   w0: PSF e^-2 radius (Gaussian PSF assumed)
%   N: number of fluorescent probes to simulate
%   k_on: rate to switch on from off-state
%   k_off: rate to switch off from on-state
%   k_p: rate to bleach
%
% OUTPUT VARIABLES
%
%   J: image series
%
% VARARGIN (NAME/VALUE PAIRS)
%
% 'kernelVarargin','kernelVar','kerVar'
% -------------------------------------------------------------------------
% {} (default) | cell-array
% -------------------------------------------------------------------------
% Varargin to pass to getImgKernel.m.
% 
% 'laserVarargin','laserVar'
% -------------------------------------------------------------------------
% {} (default) | cell-array
% -------------------------------------------------------------------------
% Varargin to pass to addLaserProfile.m. By default, no laser background
% illumination pattern is simulated.
% 
% 'model','blinkModel'
% -------------------------------------------------------------------------
% 'twoStateBleach' (default) | 'offStateBleach' | 'multiStateAllBleach'
% (see below for how to use this option)
% -------------------------------------------------------------------------
% Photo-physical/-chemical model to use. Default option is a simple on-off
% configuration and assumes equal bleaching rates from both states.
% 'offStateBleach' assumes bleaching from off-state only.
% 'multiStateAllBleach' assumes an additional on/off state. To specify the
% on/off rates of this state k_on, k_off must be supplied as 1x2 vectors
% (the second entry in each case is the corresponding on/off rate of that
% state). The value pair in this case is a cell array of form
% {'multiStateAllBleach',state_type} (see simulateMultiPhotophysics.m for
% more info). This option is not supported by Gillespie option set to 1.
%
% 'noiseType','noise'
% -------------------------------------------------------------------------
% 'emccd' (default) | 'legacy'
% -------------------------------------------------------------------------
% Noise type to include in simulation. By default, EMCCD noise is simulated
% (default values in addEMCCDNoise.m). 'legacy' is a simpler noise model
% (see addNoise.m for details).
%
% 'noiseVarargin','noiseVar'
% -------------------------------------------------------------------------
% {} (default) | cell array
% -------------------------------------------------------------------------
% Varagin to pass to noise function specified by 'noiseType' option.
%
% 'nSubFrames','subFrames','subTimeSteps'
% -------------------------------------------------------------------------
% 1 (default) | value in range [1,Inf)
% -------------------------------------------------------------------------
% This option is irrelevant if Gillespie algorithm is being used. Number of
% "sub-frames" per frame. The photo-state of each fluorophore is calculated
% at each sub-frame. A frame consists of the sum of all its sub-frames.
%
% 'offIntFrac','partOffState'
% -------------------------------------------------------------------------
% 0 (default) | value in range [0,1]
% -------------------------------------------------------------------------
% Fractional off-state intensity emission relative to on-state
%
% 'parallel', 'useParallel'
% -------------------------------------------------------------------------
% [0,0]  (default) | 1x2 vector; each entry either 0 or 1
% -------------------------------------------------------------------------
% Determines whether code should be run in parallel. First entry is for
% blinking trace generation, second entry is for image series generation
%
% 'savePath'
% -------------------------------------------------------------------------
% '' (default) | string
% -------------------------------------------------------------------------
% Save all variables in function workspace to file. Supply relative or
% absolute filename. Leave string empty (default) for no save.
%
% 'useGillespie','gillespie'
% -------------------------------------------------------------------------
% 1 (default) | 0
% -------------------------------------------------------------------------
% Set option to 1 to utilize Gillespie algorithm to simulate switch times.
% This is the most accurate way to simulate the master equation. See
% 'nSubFrames' for alternative option.
%
% -------------------------------------------------------------------------
%
% Example usage:
%
% J = blinkSim(64,2048,1.7,60,0.5,0.1,1e-4,'noisevar',...
%    {'autofluorPer',0.05},'laservar',{'laserwidth',64},'savepath',...
%    'simulation');
%
function J = blinkSim(sz,T,w0,N,k_on,k_off,k_p,varargin)

% default values for varargin options
use_parallel = [0,0];
sub_frames = 1;
blink_model = 'twoStateBleach';
use_gillespie = 1;
noise_type = 'emccd';
kernel_varargin = {};
laser_varargin = {};
noise_varargin = {};
off_int_frac = 0;
save_run = 0;
%
for i = 1:2:length(varargin)
    if any(strcmpi(varargin{i},{'parallel','useParallel'}))
        if length(varargin{i+1}) == 2
            use_parallel = varargin{i+1};
        else
            warning(['invalid option for varargin: ',varargin{i}]);
        end
    elseif any(strcmpi(varargin{i},{'nSubFrames','subFrames','subTimeSteps'}))
        if isnumeric(varargin{i+1})
            sub_frames = varargin{i+1};
        else
            warning(['invalid option for varargin: ',varargin{i}]);
        end
    elseif any(strcmpi(varargin{i},{'model','blinkModel'}))
        if iscell(varargin{i+1}) && length(varargin{i+1}) == 2 && ...
                any(strcmpi(varargin{i+1}{1},{'multiStateAllBleach'}))
            blink_model = 'multiStateAllBleach';
            state_type = varargin{i+1}{2};
        elseif any(strcmpi(varargin{i+1},{'twoStateBleach','equalBleach',...
                'eqBleach'}))
        elseif any(strcmpi(varargin{i+1},{'offStateBleach','offBleach'}))
            blink_model = 'offStateBleach';
        elseif any(strcmpi(varargin{i+1},{'asymmBleach'})) && ...
                length(k_p) == 2
            blink_model = 'asymmBleach';
        else
            warning(['Unknown option for ''',varargin{i},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{i},{'useGillespie','gillespie'}))
        if any(varargin{i+1}==[0,1])
            use_gillespie = varargin{i+1};
        else
            warning(['Unknown option for ''',varargin{i},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{i},{'noiseType','noise'}))
        if any(strcmpi(varargin{i+1},{'emccd'}))
            noise_type = 'emccd';
        elseif any(strcmpi(varargin{i+1},{'legacy'}))
            noise_type = 'legacy';
        else
            warning(['Unknown option for ''',varargin{i},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{i},{'kernelVarargin','kernelVar','kerVar'}))
        if iscell(varargin{i+1})
            kernel_varargin = varargin{i+1};
        else
            warning(['Unknown option for ''',varargin{i},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{i},{'laserVarargin','laserVar'}))
        if iscell(varargin{i+1})
            laser_varargin = varargin{i+1};
        else
            warning(['Unknown option for ''',varargin{i},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{i},{'noiseVarargin','noiseVar'}))
        if iscell(varargin{i+1})
            noise_varargin = varargin{i+1};
        else
            warning(['Unknown option for ''',varargin{i},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{i},{'offIntFrac','partOffState'}))
        if isnumeric(varargin{i+1}) && 0 <= varargin{i+1} <= 1
            off_int_frac = varargin{i+1};
        else
            warning(['Unknown option for ''',varargin{i},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{i},{'savePath'}))
        if ischar(varargin{i+1})
            save_run = 1;
            savepath = varargin{i+1};
        else
            warning(['invalid option for varargin: ',varargin{i}]);
        end
    else
        warning(['unknown varargin input ''',varargin{i},'''.'])
    end
end

% total number of frames (including sub frames)
total_T = sub_frames*T;
% interval between sub frames
sub_time = 1/sub_frames;

%% generate particle positions

position = (1 - w0) + (2*w0+sz-1).*rand(N,2);

%% photophysics

%
disp('simulating photophysics')
tic

% simulate two-state, or multi-state system
if strcmpi(blink_model,'multiStateAllBleach')
    [photo_state,obs_state,tint_obs_state] = simulateMultiPhotophysics(N,total_T,...
        sub_time,k_on,k_off,state_type,k_p,'offIntFrac',off_int_frac);
else
    if use_gillespie
        if use_parallel(1)
            [particles,tint_obs_state] = gillespiePhotophysicsPar(N,T,...
                k_on,k_off,k_p,'blinkModel',blink_model,'offIntFrac',...
                off_int_frac);
        else
            [init_photo_state,state_changes,tint_obs_state,n_cycles] = ...
                gillespiePhotophysics(N,T,k_on,k_off,k_p,'blinkModel',...
                blink_model,'offIntFrac',off_int_frac);
        end
    else
        [photo_state,obs_state,tint_obs_state] = simulatePhotophysics(N,total_T,sub_time,k_on,k_off,k_p,...
            'blinkModel',blink_model,'offIntFrac',off_int_frac);
    end
end

toc
%

%% image series creation

% array to store noiseless image series
J_sig = zeros(sz,sz,T);

%
disp('generating image series')
tic

img_kernel_stat = zeros(sz,sz,N);
% get corresponding image kernel
if w0 == 0
    for i = 1:N
        round_y = mod(round(position(i,1)),sz)+1;
        round_x = mod(round(position(i,2)),sz)+1;
        img_kernel_stat(round_y,round_x,i) = 1;
    end
else
    for i = 1:N
        img_kernel_stat(:,:,i) = getImgKernel(J_sig,position(i,:),w0,...
            kernel_varargin{:});
    end
end

if use_parallel(2)
    parfor t = 1:T
        % time-integrated image at time t
        obs_kernel_stat = tint_obs_state(t,1,:).*sub_time.*img_kernel_stat;
        J_sig(:,:,t) = sum(obs_kernel_stat,3);
    end
    delete(gcp)
else
    for t = 1:T
        % time-integrated image at time t
        obs_kernel_stat = tint_obs_state(t,1,:).*sub_time.*img_kernel_stat;
        J_sig(:,:,t) = sum(obs_kernel_stat,3);
    end
end
toc
%

%
disp('generating noise')
tic

if strcmp(noise_type,'emccd')
    J = addEMCCDNoise(J_sig,laser_varargin,noise_varargin{:});
else
    [J,shot_noise_ser,wg_noise_ser] = addNoise(J_sig,noise_varargin{:});
end

toc
%

% prevent double saving of noiseless movies
if all(J(:) == J_sig(:))
    clear J_sig
end

%% save simulation info

if save_run
    [save_dir,~,save_ext] = fileparts(savepath);
    % prepend current directory path if savepath empty
    if isempty(save_dir)
        savepath = [cd,filesep,savepath];
    end
    % append default extension '.mat' to savepath if empty
    if isempty(save_ext)
        savepath = [savepath,'.mat'];
    end
    %
    if exist(savepath,'file')
        warning('file already exists; not overwritten')
    else
        %     [filepath,filename] = fileparts(savepath);
        save(savepath,'-v7.3')
    end
end