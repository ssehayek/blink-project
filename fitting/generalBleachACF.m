% ACF of blinking and independent bleaching rates from on/off states (with
% time-integration). See demo for usage when fitting.
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
function [out] = generalBleachACF(params,tau_vec,n_frames,varargin)

% initialize no variables as symbolic
all_vars = {'a','b','k1','k2','T','tau'};
all_vals = {params(1),params(2),params(3),params(4),n_frames,tau_vec};
sym_vars = {''};
s = setSymVars(all_vars,all_vals,sym_vars);
%

for ii = 1:length(varargin)
    if ~iscell(varargin{ii})
        % choose symbolic variables
        if any(strcmpi(varargin{ii},{'symVars','symsVars'}))
            if iscell(varargin{ii+1})
                sym_vars = varargin{ii+1};
                s = setSymVars(all_vars,all_vals,sym_vars);
            else
                warning(['Unknown option for ''',varargin{ii},...
                    ''', using default options.'])
            end
        end
    end
end

K = s.a+s.b;
p0 = s.a/K;
Kp = K + s.k1+s.k2;
d = sqrt(Kp^2-4*(s.a*s.k1+(s.b+s.k1).*s.k2));

L1 = (Kp+d)/2;
L2 = (Kp-d)/2;

v_11 = -(Kp+d-2*s.k2)/(2*(s.k1-s.k2));
v_12 = (Kp+d-2*s.k1)/(2*(s.k1-s.k2));
v_21 = -(Kp-d-2*s.k2)/(2*(s.k1-s.k2));
v_22 = (Kp-d-2*s.k1)/(2*(s.k1-s.k2));

a_1 = v_22/(v_11*v_22-v_12*v_21);
a_2 = v_12/(v_12*v_21-v_11*v_22);

b_1 = -a_1/v_22*v_21;
b_2 = -a_2/v_12*v_11;

c_1 = v_11*(a_1*p0+b_1*(1-p0));
c_2 = v_21*(a_2*p0+b_2*(1-p0));

D1 = c_1*a_1*v_11;
D2 = c_1*a_2*v_21;
D3 = c_2*a_1*v_11;
D4 = c_2*a_2*v_21;


Phi_tau_num = ((-1)+exp(1).^L1).^(-1).*((-1)+exp(1).^L2).^(-1).*L1.^(-1).*(L1+( ...
    -1).*L2).^(-1).*L2.^(-1).*((-1).*D2.*exp(1).^((-1).*L2.*(1+s.tau)).*(( ...
    -1)+exp(1).^L2).^2.*((-1).*exp(1).^L1+exp(1).^L2).*(1+(-1).*exp(1) ...
    .^(L1+(-1).*L1.*(s.T+(-1).*s.tau))).*L1+((-1)+exp(1).^L1).*(D4.*exp(1) ...
    .^((-1).*L2.*s.tau).*((-1)+exp(1).^L2).*(1+(-1).*exp(1).^((-1).*L2.*( ...
    s.T+(-1).*s.tau))).*L1.*(L1+(-1).*L2)+(D3.*exp(1).^((-1).*L1.*(1+s.tau)).*(( ...
    -1)+exp(1).^L1).*(exp(1).^L1+(-1).*exp(1).^L2).*(1+(-1).*exp(1).^( ...
    (-1).*L2.*(s.T+(-1).*s.tau)))+D1.*exp(1).^((-1).*L1.*s.tau).*((-1)+exp(1) ...
    .^L2).*(1+(-1).*exp(1).^((-1).*L1.*(s.T+(-1).*s.tau))).*(L1+(-1).*L2)).* ...
    L2)).*(s.T+(-1).*s.tau).^(-1);
Phi_tau_denom = ((-1)+exp(1).^L1).^(-1).*((-1)+exp(1).^L2).^(-1).*L1.^(-1).*(L1+( ...
    -1).*L2).^(-1).*L2.^(-1).*((-1).*D2.*exp(1).^((-1).*L2.*2).*(( ...
    -1)+exp(1).^L2).^2.*((-1).*exp(1).^L1+exp(1).^L2).*(1+(-1).*exp(1) ...
    .^(L1+(-1).*L1.*(s.T-1))).*L1+((-1)+exp(1).^L1).*(D4.*exp(1) ...
    .^((-1).*L2).*((-1)+exp(1).^L2).*(1+(-1).*exp(1).^((-1).*L2.*( ...
    s.T-1))).*L1.*(L1+(-1).*L2)+(D3.*exp(1).^((-1).*L1.*2).*(( ...
    -1)+exp(1).^L1).*(exp(1).^L1+(-1).*exp(1).^L2).*(1+(-1).*exp(1).^( ...
    (-1).*L2.*(s.T-1)))+D1.*exp(1).^((-1).*L1).*((-1)+exp(1) ...
    .^L2).*(1+(-1).*exp(1).^((-1).*L1.*(s.T-1))).*(L1+(-1).*L2)).* ...
    L2)).*(s.T-1).^(-1);
Phi_tau = Phi_tau_num./Phi_tau_denom;

out = double(Phi_tau);