function Td = discretize_opvar(T,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Td] = discretize_opvar(T,N) takes in a 4-PI operator T that acts on
% functions in the interval [a,b] and changes it to a matrix Td that acts
% on polynomials of degree decomposed using chebyshev polynomial basis upto
% degree N.
%
% INPUT:
%
% T: a 4-PI operator
% N: max degree of chebyshev polynomial basis, default 2
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - discretize_opvar
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
% Initial coding MMP, SS  - 1_29_2021


if nargin==1
    N = 2;
elseif nargin>2
    error("Incorrect number of inputs; only 2 inputs are allowed");
end
if ~isa(T,'opvar')
    error("First argument must be a opvar class object");
end
if mod(N,1)~=0
    error("Second argument must be an integer");
end

I_init = T.I;
a = I_init(1);
b = I_init(2);

Tn = transl(T,[-1,1]);

s = Tn.var1;
theta = Tn.var2;

dim = Tn.dim;

%get dim of input vectors
m = dim(1,2); n = dim(2,2);

%
u_f = polynomial(zeros(m+n,1)); 

%generate pvar coeffs
if m>0
for i=1:m
    eval(['pvar cf_' num2str(i)]);
    u_f(i,1) = eval(['cf_' num2str(i)]);
end
for i=1:n
    for j =0:N
        eval(['pvar cf_' num2str(m+j+(i-1)*N)]);
        u_f(m+i,1) = u_f(i,1) + eval(['cf_' num2str(m+j+(i-1)*N)])*cheb_poly_2(s,j);
    end
end
elseif m==0
for i=1:n
    for j =0:N
        eval(['pvar cf_' num2str(j+(i-1)*N)]);
        u_f(i,1) = u_f(i,1) + eval(['cf_' num2str(j+(i-1)*N)])*cheb_poly_2(s,j);
    end
end
end
% v_f =  Tu_f

v_f = [Tn.P*u_f(1:m)+int(Tn.Q1*u_f(m+1:m+n),s,-1,1); 
       Tn.Q2*u_f(1:m)+Tn.R.R0*u_f(m+1:m+n)+int(Tn.R.R1*u_f(m+1:m+n),theta,-1,s)+int(Tn.R.R2*u_f(m+1:m+n),theta,s,1)];

end

