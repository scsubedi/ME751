function varargout = getdp1(gcon_inp,varargin)
% function that calculates phi, nu,gamma, dphi_dr, dphi_dp for GCon DP1


%%%%%%%%%%%% enter vararg conditions here %%%%%%%%%%%%%%
% extracting input variables from structure
a_1 = gcon_inp.a_1;       % given vector in local coordinates, i.e. a1bar
p1 = gcon_inp.p1;
p1dot = gcon_inp.p1dot;
a_2 = gcon_inp.a_2;
p2 = gcon_inp.p2;
p2dot = gcon_inp.p2dot;
ft = gcon_inp.ft;
dft = gcon_inp.dft;
ddft = gcon_inp.ddft;
i = gcon_inp.i;
j = gcon_inp.j;

A1 = getA(p1);
    a1 = A1*a_1;        % finding vector in global coordinates
A2 = getA(p2);
    a2 = A2*a_2;

% calculating all the B matrices required
B_p2 = getB(p2,a_2);
B_p2dot = getB(p2dot,a_2);   
B_p1 = getB(p1,a_1);
B_p1dot = getB(p1dot,a_1);
    
phi = a_1.'*A1.'*A2*a_2 - ft;

nu = dft;

gamma = -a1.'*B_p2dot*p2dot - a2.'*B_p1dot*p1dot - ...
    2*(B_p1*p1dot).'*(B_p2*p2dot) + ddft;

ddr1 = zeros(1,3);               % DP1 has no r terms
ddr2 = zeros(1,3);
ddp1 = a2.'*B_p1;
ddp2 = a1.'*B_p2;

% checking if either bodies are ground
if i == 0
    dphi_dr = ddr2;
    dphi_dp = ddp2;
elseif j == 0
    dphi_dr = ddr1;
    dphi_dp = ddp1;
else
    dphi_dr = [ddr1,ddr2];
    dphi_dp = [ddp1,ddp2];
end

arg = 1;
if nargin > 1
    if any(strcmp('phi',varargin{1}))
        varargout{arg} = phi;
        arg = arg+1;
    end
    if any(strcmp('nu',varargin{1}))
        varargout{arg} = nu;
        arg = arg+1;
    end
    if any(strcmp('gamma',varargin{1}))
        varargout{arg} = gamma;
        arg = arg+1;
    end
    if any(strcmp('dphi_dr',varargin{1}))
        varargout{arg} = dphi_dr;
        arg = arg+1;
    end
    if any(strcmp('dphi_dp',varargin{1}))
        varargout{arg} = dphi_dp;
    end
else
    varargout = {phi,nu,gamma,dphi_dr,dphi_dp};
end
    


