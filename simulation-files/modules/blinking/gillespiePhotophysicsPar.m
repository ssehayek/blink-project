function [particles,tint_obs_state] = gillespiePhotophysicsPar(N,T,k_on,k_off,k_p,varargin)

off_int_frac = 0;
blink_model = 'twoStateBleach';
for ii = 1:2:length(varargin)
    if any(strcmpi(varargin{ii},{'offIntFrac','partOffState'}))
        if isnumeric(varargin{ii+1}) && 0 <= varargin{ii+1} && ...
                varargin{ii+1} <= 1
            off_int_frac = varargin{ii+1};
        end
    elseif any(strcmpi(varargin{ii},{'model','blinkModel'}))
        if any(strcmpi(varargin{ii+1},{'twoStateBleach','equalBleach',...
                'eqBleach'}))
            % blinking + bleaching from both states (default)
        elseif any(strcmpi(varargin{ii+1},{'offStateBleach','offBleach'}))
            % blinking + bleaching from off-state
            blink_model = 'offStateBleach';
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    end
end

%% initial states

% sum of rates
K = k_on + k_off;
% expected time between reactions (ignoring bleaching)
mean_dt_est = (2*k_on*k_off/K)^(-1);
% expected number of reactions (ignoring bleaching)
mean_flips_est = round(T/mean_dt_est);

% photo-state initialization
%
init_photo_state = 1 + binornd(1,k_off/K,[N,1]);
%
init_obs_state = (init_photo_state == 1);
% convert logical to double
init_obs_state = double(init_obs_state);
% assign off-state "intensity" to off_int_frac
init_obs_state(init_obs_state == 0) = off_int_frac;
%

%% state update using Gillespie algorithm

switch_time_temp = cell(1,N);
photo_state_temp = cell(1,N);
obs_state_temp = cell(1,N);
for n = 1:N
    % initialize according to mean number of switches
    %
    switch_time_temp{n} = nan(1,mean_flips_est);
    photo_state_temp{n} = nan(1,mean_flips_est);
    obs_state_temp{n} = nan(1,mean_flips_est);
    % initial values
    %
    switch_time_temp{n}(1) = 0;
    photo_state_temp{n}(1) = init_photo_state(n);
    obs_state_temp{n}(1) = init_obs_state(n);
end

for n = 1:N
    % intialize switch counter for particle n
    n_s = 0;
    % intialize time for particle n
    t = 0;
    % sliced 
    switch_time_temp_n = switch_time_temp{n};
    photo_state_temp_n = photo_state_temp{n};
    obs_state_temp_n = obs_state_temp{n};
    
    % loop breaks if t > T, or particle bleaches
    while true
        % current photo-state of particle n
        s = photo_state_temp_n(n_s+1);
        % array of propensities
        switch blink_model
            case 'twoStateBleach'
                a_j = [(s==2)*k_on,(s==1)*k_off,k_p];
            otherwise
                a_j = [(s==2)*k_on,(s==1)*k_off,(s==2)*k_p];
        end
        % cumulative sum of propensities
        S_a_j = cumsum(a_j);
        % sum of propensities
        a_0 = S_a_j(end);
        
        %
        r_1 = rand();
        r_2 = rand();
        
        % time to next state switch
        dt = -1/a_0*log(r_1);
        % next state
        j = find(S_a_j >= r_2*a_0,1,'first');
        
        % time update
        t = t + dt;
        
        if t <= T
            % increase switch counter
            n_s = n_s + 1;
        else
            % break if t > T
            break
        end
        
        % update particles struct
        switch_time_temp_n(n_s+1) = t;
        photo_state_temp_n(n_s+1) = j;
        obs_state_temp_n(n_s+1) = (j==1)+(j==2)*off_int_frac;
        
        if j == 3
            % break if particle bleaches
            break
        end
    end
    % find first index that contains nan (if it exists)
    nan_ind = find(isnan(switch_time_temp_n),1,'first');
    % delete all indices after nan_ind
    switch_time_temp_n(nan_ind:end) = [];
    photo_state_temp_n(nan_ind:end) = [];
    obs_state_temp_n(nan_ind:end) = [];
        
    switch_time_temp{n} = switch_time_temp_n;
    photo_state_temp{n} = photo_state_temp_n;
    obs_state_temp{n} = obs_state_temp_n;
end

%% time-integrated observed states

% time-integrated "intensity" of each particle as a function of frame.
% a complete array listing all states after each switch is impractical as
% it can exceed size limitations and contains redundant information.
% useful for later image series construction.
tint_obs_state = zeros(N,T);
parfor n = 1:N
    % state switch times for nth molecule
    switch_times_n = switch_time_temp{n};
    % floor of switch times; useful for classifying which frame each
    % switch belongs to
    switch_times_flr_n = floor(switch_times_n);
    
    % observed states of nth particle
    obs_states_n = obs_state_temp{n};
    
    % temporary variable for storing observed state at the beginning of
    % next frame
    next_start_state = obs_states_n(1);
    
    %
    tint_obs_state_n = tint_obs_state(n,:);
    % frame loop
    for t = 0:T-1
        % find indices of switches which occur in frame t
        switch_inds_t = (switch_times_flr_n == t);
        
        % observed states at each switch in frame t (including state at
        % start of frame)
        curr_obs_state_t = [next_start_state,obs_states_n(switch_inds_t)];
        
        % time at which switches occur in frame t (including start and end
        % times of frame)
        switch_times_t = [t,switch_times_n(switch_inds_t),t+1];
        % time between consecutive switches
        d_switch_time_t = diff(switch_times_t);
        
        % store time-integrated observed state
        tint_obs_state_n(t+1) = dot(curr_obs_state_t,d_switch_time_t);
        
        % update state at the beginning of next frame
        next_start_state = curr_obs_state_t(end);
    end
    tint_obs_state(n,:) = tint_obs_state_n;
end
% permute dimensions for later multiplication with img_kernel_stat
tint_obs_state = permute(tint_obs_state,[2,3,1]);

% particles is the output array which stores info about states and switch
% times
particles.switch_time = switch_time_temp;
particles.photo_state = photo_state_temp;
particles.obs_state = obs_state_temp;

delete(gcp)