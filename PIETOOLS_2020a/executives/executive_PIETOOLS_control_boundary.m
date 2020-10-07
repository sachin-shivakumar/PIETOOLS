%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script executes a synthesis code for H-infty optimal controller 
% design (control at the boundary, not in domain) for a 4-PIE System defined
% by the 9 4-PI operator representation
% Eop \dot x = Aop x(t)+B1op w(t)+B2op u(t)
% z(t)=C1op x(t)+D11 w(t)+D12 u(t)
% y(t)=C2op x(t)+D21 w(t)+D22 u(t)
% u(t)=Kop x(t)
%
% Where for now, we assume $z(t)\in \R^{nz}$, $x(t) \in \R^{nx1} \times
% L_2^{nx2}$, w(t)\in \R^{nw}, u(t)\in \R^{nu}, y(t) \in \R^{ny}
% The domain is s \in [a,b]
%
% The following must be defined externally:
% Eop,Aop,B1op,C1op,B2op,C2op,D11,D12,D21,D22,nu,ny,nw,nz,nx1,nx2,n_order1,n_order2,a,b
% s and th must be pvars
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PIETools - executive_PIETOOLS_Hinf_control
%
% Copyright (C)2019  M. Peet, S. Shivakumar
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, SS  - 7_26_2019
%

if ~exist('override1')
    override1=1; % default to no Psatz in the LF
end

if ~exist('override2')
    override2=0; % default to Psatz in the Derivative
end

% setup an SOS program

varlist = [s; theta];
prog = sosprogram(varlist);

% domain of spatial variables
% p-compact form of the interval
%g1 = (X(2)-s)*(s-X(1));

% hinf norm variable which needs to be minimized
%prog = sosdecvar(prog, gamma); %this sets gamma as decision var
% prog = sossetobj(prog, gamma); %this minimizes gamma, comment for feasibility test
options1.sep = 1; options12.sep = 1;
% setup the variables for lyapunov function candidate
disp('Parameterizing Positive function');
[prog, P1op] = poslpivar(prog, [nx1 ,nx2],X,dd1,options1);

if override1~=1
    [prog, P2op] = poslpivar(prog, [nx1 ,nx2],X,dd12,options12);
    Pop=P1op+P2op;
else
    Pop=P1op;
end

% enforce strict positivity on the operator
Pop.P = Pop.P+eppos*eye(nx1);
Pop.R.R0 = Pop.R.R0+eppos2*eye(nx2);  

[prog,Zop] = lpivar(prog,[nu nx1;0 nx2],X,ddZ);


% Initialize zero operators
opvar ZZ;
ZZ.dim = [nx1 nx1; nx2 nx2]; ZZ.I = [a b];

%Assemble the big operator
Dop =   [-Pop                        ZZ         Zop'*B2op';
         ZZ                          -Pop       Zop'*TB2op';
         B2op*Zop                     TB2op*Zop     Top*(B2op*Zop+Aop*Pop)'+(Aop*Pop+B2op*Zop)*Top'+TB2op*Zop*Aop'+Aop*Zop'*TB2op']; 

disp('Parameterize the derivative inequality');
disp('Enforcing the Negativity Constraint');
opts.pure=0;

if sosineq_on
    disp('Using lpi_ineq');
    prog = lpi_ineq(prog,-Dop,opts);
else
    disp('Using an Equality constraint');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% The old way, where you have to specify everything
    [prog, De1op] = poslpivar(prog, [nw+nz+nx1, nx2],X,dd2,options2);
    
    if override2~=1
        [prog, De2op] = poslpivar(prog,[nw+nz+nx1, nx2],X, dd3,options3);
        Deop=De1op+De2op;
    else
        Deop=De1op;
    end
    % derivative negativity
    % constraints
    prog = lpi_eq(prog,Deop+Dop); %Dop=-Deop
end

%solving the sos program
prog = sossolve(prog,sos_opts); 
