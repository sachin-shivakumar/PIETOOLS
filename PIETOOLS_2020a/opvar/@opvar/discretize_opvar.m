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
if isinteger(N)~=1
    error("Second argument must be an integer");
end

I_init = T.I;
a = I_init(1);
b = I_init(2);

c = I(1);
d = I(2);

pvar s theta;

end