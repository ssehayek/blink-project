function [init_photo_state,state_changes,tint_obs_state,n_cycles] = gillespiePhotophysics(N,T,k_on,k_off,k_p,varargin)

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
        elseif any(strcmpi(varargin{ii+1},{'asymmBleach'})) && ...
                length(k_p) == 2
            blink_model = 'asymmBleach';
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
mean_dt_est = (2*N*k_on*k_off/K)^(-1);
% expected number of reactions (ignoring bleaching)
mean_flips_est = round(T/mean_dt_est);

% p_bleach = expcdf(T,1/k_p);
% mean_bleach_est = round(N*p_bleach);

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

% absolute time
time = 0;
% state_changes stores info about state changes, and has size based on
% estimated number of switches
%
% first row: indicates particle index which undergoes switch
% second row: time of switch
% third row: new state at this time
% fourth row: observed state at this time
state_changes = nan(4,mean_flips_est);

% temporary array for storing current particle photo-states
photo_state = init_photo_state;

% initial number of on-state particles
N_on = nnz(init_photo_state==1);
% initial number of off-state particles
N_off = N - N_on;

% counter for number of switches
n = 0;
% array for storing number of state cycles for each molecule (excluding
% bleaching)
% each state switch counts for 0.5 of a cycle
n_cycles = zeros(N,1);
while time <= T-1 && N_on + N_off > 0
    % array of propensities
    switch blink_model
        case 'twoStateBleach'
            a_j = [N_off*k_on,N_on*k_off,(N_on+N_off)*k_p];
        case 'offStateBleach'
            a_j = [N_off*k_on,N_on*k_off,N_off*k_p];
        otherwise
            a_j = [N_off*k_on,N_on*k_off,N_on*k_p(1)+N_off*k_p(2)];
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
    time = time + dt;
    % increase switch counter
    n = n + 1;
    switch j
        case 1
            % get all particle indices which are in off-state
            prev_state_inds = find(photo_state==2);
            % index of random off particle to switch on
            switch_mol_ind = prev_state_inds(randsample(length(prev_state_inds),1));
            
            % update particle numbers
            N_on = N_on + 1;
            N_off = N_off - 1;
            
            % increase blink cycle counter
            d_cycle = k_off/(k_on+k_off);
            n_cycles(switch_mol_ind) = n_cycles(switch_mol_ind) + d_cycle;
        case 2
            % get all particle indices which are in on-state
            prev_state_inds = find(photo_state==1);
            % index of random on particle to switch off
            switch_mol_ind = prev_state_inds(randsample(length(prev_state_inds),1));
            
            % update particle numbers
            N_on = N_on - 1;
            N_off = N_off + 1;
            
            % increase blink cycle counter
            d_cycle = k_on/(k_on+k_off);
            n_cycles(switch_mol_ind) = n_cycles(switch_mol_ind) + d_cycle;
        case 3
            switch blink_model
                case 'twoStateBleach'
                    % get all particle indices which are either on or off
                    prev_state_inds = find(photo_state~=3);
                    % random particle index to bleach
                    switch_mol_ind = prev_state_inds(randsample(length(prev_state_inds),1));
                    
                    % update particle numbers depending on the previously occupied
                    % state
                    prev_state = photo_state(switch_mol_ind);
                    N_on = N_on - (2-prev_state);
                    N_off = N_off - (prev_state-1);
                case 'offStateBleach'
                    % get all particle indices which are in off-state
                    prev_state_inds = find(photo_state==2);
                    % index of random off particle to bleach
                    switch_mol_ind = prev_state_inds(randsample(length(prev_state_inds),1));
                    
                    % update particle numbers
                    N_off = N_off - 1;
                otherwise
                    % first randomly draw which state will bleach
                    % 
                    % probability that on-state bleaches
                    p_on_bleach = N_on*k_p(1)/a_j(end);
                    % state to bleach 
                    state_bleach = 2-binornd(1,p_on_bleach);
                    %
                    % get all particle indices which are in state_bleach
                    prev_state_inds = find(photo_state==state_bleach);
                    % index of random particle in the state state_bleach to bleach
                    switch_mol_ind = prev_state_inds(randsample(length(prev_state_inds),1));
                    %
                    % update particle numbers depending on the previously
                    % occupied state
                    N_on = N_on - (2-state_bleach);
                    N_off = N_off - (state_bleach-1);
            end
    end
    % update state_changes based on current switch
    state_changes(1,n) = switch_mol_ind;
    state_changes(2,n) = time;
    state_changes(3,n) = j;
    state_changes(4,n) = (j==1)+(j==2)*off_int_frac;
    
    photo_state(switch_mol_ind) = j;
end
% find first column that contains nan (if it exists)
[~,nan_col]=find(isnan(state_changes),1,'first');
% delete all columns after this column
state_changes(:,nan_col:end) = [];

%% time-integrated observed states

% time-integrated "intensity" of each particle as a function of frame.
% a complete array listing all states after each switch is impractical as
% it can exceed size limitations and contains redundant information.
% useful for later image series construction.
tint_obs_state = zeros(N,T);
for n = 1:N
    % state switch times for nth molecule
    switch_times = state_changes(2,state_changes(1,:) == n);
    % floor of switch times; useful for classifying which frame each
    % switch belongs to
    switch_times_flr = floor(switch_times);
    
    % observed states of nth particle
    curr_obs_state = state_changes(4,state_changes(1,:) == n);
    
    % temporary variable for storing state at the beginning of next frame
    next_start_state = init_obs_state(n);
    % frame loop
    for t = 0:T-1
        % find indices of switches which occur in frame t
        switch_inds_t = (switch_times_flr == t);
        
        % observed states at each switch in frame t (including state at
        % start of frame)
        curr_obs_state_t = [next_start_state,...
            curr_obs_state(switch_inds_t)];
        
        % time at which switches occur in frame t (including start and end
        % times of frame)
        switch_times_t = [t,switch_times(switch_inds_t),t+1];
        % time between consecutive switches
        d_switch_times_t = diff(switch_times_t);
        
        % store time-integrated observed state
        tint_obs_state(n,t+1) = dot(curr_obs_state_t,d_switch_times_t);
        
        % update state at the beginning of next frame
        next_start_state = curr_obs_state_t(end);
    end
end
% permute dimensions for later multiplication with img_kernel_stat
tint_obs_state = permute(tint_obs_state,[2,3,1]);