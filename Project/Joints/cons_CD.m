function varargout = cons_CD(gcon_inp,varargin)
% function that calculates phi, nu,gamma, dphi_dr, dphi_dp for GCon CD


%%%%%%%%%%%% enter vararg conditions here %%%%%%%%%%%%%%
% extracting input variables from structure
s1_p = gcon_inp.s1_p;       % given vector in local coordinates, i.e. a1bar
p1 = gcon_inp.p1;
p1dot = gcon_inp.p1dot;
s2_q = gcon_inp.s2_q;
p2 = gcon_inp.p2;
p2dot = gcon_inp.p2dot;
c = gcon_inp.c;
r1 = gcon_inp.r1;
r2 = gcon_inp.r2;
ft = gcon_inp.ft;
dft = gcon_inp.dft;
ddft = gcon_inp.ddft;
i = gcon_inp.i;
j = gcon_inp.j;

A1 = Amatrix(gcon_inp.p1);
s1p = A1*gcon_inp.s1_p;        % finding vector in global coordinates
A2 = Amatrix(gcon_inp.p2);
s2q = A2*gcon_inp.s2_q;

% calculating all the B matrices required
B_p2 = Bmatrix(gcon_inp.p2,gcon_inp.s2_q);
B_p2dot = Bmatrix(gcon_inp.p2dot,gcon_inp.s2_q);
B_p1 = Bmatrix(gcon_inp.p1,gcon_inp.s1_p);
B_p1dot = Bmatrix(gcon_inp.p1dot,gcon_inp.s1_p);

phi = gcon_inp.c'*(gcon_inp.r2 + gcon_inp.s2_q - gcon_inp.r1 - gcon_inp.s1_p) - gcon_inp.ft;

nu = gcon_inp.dft;

gamma = gcon_inp.c.'*B_p1dot*gcon_inp.p1dot - gcon_inp.c.'*B_p2dot*gcon_inp.p2dot + gcon_inp.ddft;

ddr1 = -gcon_inp.c.'*eye(3);               % DP1 has no r terms
ddr2 = gcon_inp.c.'*eye(3);
ddp1 = -gcon_inp.c.'*B_p1;
ddp2 = gcon_inp.c.'*B_p2;

% checking if either bodies are ground
if gcon_inp.i == 0
    dphi_dr = ddr2;
    dphi_dp = ddp2;
elseif gcon_inp.j == 0
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