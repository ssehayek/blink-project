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
init_photo_state = 1 + binornd(1,k_off/K,[N,1]);

%% state update using Gillespie algorithm

% absolute time
t = 0;
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
while t <= T-1 && N_on + N_off > 0
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
    t = t + dt;
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
                    prev_state = 2-binornd(1,p_on_bleach);
                    %
                    % get all particle indices which are in state_bleach
                    prev_state_inds = find(photo_state==prev_state);
                    % index of random particle in the state state_bleach to bleach
                    switch_mol_ind = prev_state_inds(randsample(length(prev_state_inds),1));
                    %
                    % update particle numbers depending on the previously
                    % occupied state
                    N_on = N_on - (2-prev_state);
                    N_off = N_off - (prev_state-1);
            end
    end
    % update state_changes based on current switch
    state_changes(1,n) = switch_mol_ind;
    state_changes(2,n) = t;
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
% possible values of observed state for ease of computation
obs_state_vals = [1,off_int_frac];
init_obs_state = (init_photo_state == 1) + off_int_frac.*(init_photo_state == 2);
for n = 1:N
    % state switches in nth dye
    switch_times_n = state_changes(2,state_changes(1,:)==n);
    % index when nth fluorophore bleaches
    t_bleach_ind = find(state_changes(3,state_changes(1,:)==n)==3,1,'first');
    % time when nth fluorophore bleaches
    t_bleach_n = switch_times_n(t_bleach_ind);
    %
    % number of frames until nth fluorophore bleaches; T if no bleach
    T_eff = min([ceil(t_bleach_n),T]);
    
    frame_vec = 0:T_eff;
    % times at which switches occur for particle n (including start and end
    % frame times))
    switch_times_n = sort([frame_vec,state_changes(2,state_changes(1,:)==n)]);
    % find frame indices in sorted switch_times_n
    frame_inds = find(switch_times_n == ceil(switch_times_n));
    % temporary variable for storing observed/photo-state at the beginning
    % of next frame
    next_photo_state = init_photo_state(n);
    next_obs_state = init_obs_state(n);
    % frame loop
    for t = 1:T_eff-1
        % get switch times for nth fluorophore in frame t
        switch_times_n_t = switch_times_n(frame_inds(t):frame_inds(t+1));
        % number of switches in frame t
        n_switch = length(switch_times_n(frame_inds(t):frame_inds(t+1)))-2;
        % times between consecutive switches
        d_switch_times_t = diff(switch_times_n_t);
        if n_switch == 0 && next_obs_state == 0
            tint_obs_state(n,t) = 0;
        elseif off_int_frac == 0 && next_obs_state == 1
            tint_obs_state(n,t) = sum(d_switch_times_t(1:2:end));
        elseif off_int_frac == 0 && next_obs_state == 0
            tint_obs_state(n,t) = sum(d_switch_times_t(2:2:end));
        elseif next_obs_state == 1
            % sum times of odd switch intervals
            dt_odd = sum(d_switch_times_t(1:2:end));
            % store time-integrated observed state in last frame (faster
            % than using built-in dot function)
            tint_obs_state(n,t) = dt_odd + (dt_odd-1)*off_int_frac;
        elseif next_obs_state == off_int_frac
            dt_odd = sum(d_switch_times_t(1:2:end));
            tint_obs_state(n,t) = dt_odd*off_int_frac + (dt_odd-1);
        end
        
        if mod(n_switch,2) == 1
            % if n_switch is odd, switch next_photo_state and
            % next_obs_state
            next_photo_state = 3-next_photo_state;
            next_obs_state = obs_state_vals(next_photo_state);
        end
    end
    % last frame (sacrifice code efficiency for compactness)
    %
    % get switch times for nth fluorophore in last frame
    switch_times_n_t = switch_times_n(frame_inds(T_eff):frame_inds(T_eff+1));
    % times between consecutive switches in last frame
    d_switch_times_t = diff(switch_times_n_t);
    %
    if isempty(t_bleach_n)
        % fluorophore n does not bleach
        %
        dt_odd = sum(d_switch_times_t(1:2:end));
        % store time-integrated observed state in last frame
        tint_obs_state(n,T_eff) = dt_odd*next_obs_state+(1-dt_odd)*...
            obs_state_vals(obs_state_vals ~= next_obs_state);
    else
        % fluorophore n bleaches
        %
        dt_odd = sum(d_switch_times_t(1:2:end-1));
        dt_even = sum(d_switch_times_t(2:2:end-1));
        %
        tint_obs_state(n,T_eff) = dt_odd*next_obs_state+dt_even*...
            obs_state_vals(obs_state_vals ~= next_obs_state);
    end
end
% permute dimensions for later multiplication with img_kernel_stat
tint_obs_state = permute(tint_obs_state,[2,3,1]);