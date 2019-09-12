% ACF of blinking and equal bleaching rates from on/off states (with
% time-integration). See demo for usage when fitting. If bleaching is not
% significant, one can avoid fitting for bleaching by specifying the fit
% function as: 
% fit_fun = @(p,tau) equalBleachACF([p(1:2),eps],tau,length(phi_tau));
%
% INPUT VARIABLES
%
% params: fit parameters (kon,koff,kp); if fitting, this is left as an
%         input argument for the function handle
% tau_vec: time-lags considered in ACF; if fitting, this is left as an
%          input argument for the function handle
% n_frames: number of frames considered
% 
% VARARGIN
% 
% 'symVars', 'symsVars' (not recommended)
% -------------------------------------------------------------------------
% {''} (default) | any subset of {'a','b','kp','T','tau'}
% -------------------------------------------------------------------------
% Option to use Symbolic Toolbox when calculating ACF. The variables
% specified in this option will be treated as symbolic before evaluation.
%
function [out] = equalBleachACF(params,tau_vec,n_frames,...
    varargin)

% ordered variable names
all_vars = {'a','b','kp','T','tau'};
all_vals = {params(1),params(2),params(3),n_frames,tau_vec};
%

sym_vars = {''}; % initialize no variables as symbolic
err_bool = 0;
for ii = 1:length(varargin)
    if ~iscell(varargin{ii})
        if any(strcmpi(varargin{ii},{'symVars','symsVars'}))
            % choose symbolic variables
            if iscell(varargin{ii+1}) && all(ischar([varargin{ii+1}{:}]))
                sym_vars = varargin{ii+1};
            else
                warning(['Unknown option for ''',varargin{ii},...
                    ''', using default options.'])
            end
        elseif any(strcmpi(varargin{ii},{'error','err','residual','res'}))
            err_bool = 1;
            ydata = varargin{ii+1};
        end
    end
end
s = setSymVars(all_vars,all_vals,sym_vars);

K = s.a + s.b;
K_p = K + s.kp;
Phi_tau_num = 1./(s.T-s.tau).*(exp(-s.kp.*s.tau)-exp(-s.kp.*s.T)).*...
    (s.b.*s.kp.*(1-exp(-K)).*(1-exp(-K_p)).*exp(-K.*(s.tau-1))+...
    s.a.*K.*K_p.*(1-exp(-s.kp)));
Phi_tau_denom = 1./(s.T-1).*(exp(-s.kp)-exp(-s.kp.*s.T)).*...
    (s.b.*s.kp.*(1-exp(-K)).*(1-exp(-K_p))+...
    s.a.*K.*K_p.*(1-exp(-s.kp)));
Phi_tau = Phi_tau_num./Phi_tau_denom;

switch err_bool
    case 1
        ydata = reshape(ydata,size(Phi_tau));
        err = norm(Phi_tau-ydata)^2;
        out = double(err);
    otherwise
        out = double(Phi_tau);
end