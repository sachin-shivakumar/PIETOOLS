clc; clear;
% Random ODE-PDE generator
% pvar s theta;
% ni = {1,2,3,5}; N = length(ni)-1;
% nz = 1; ny = 2; nw =2; nu=1; nx = 2;
% a = 0; b = 1;
% 
% np = 0; nrL1 = 0;
% for i =0:N
% np=np+ni{i+1};
% nrL1 = nrL1+i*ni{i+1};
% end
% nrL1 = 2*nrL1;
% 
% A = rand(nx);
% B11 = rand(nx,nw); B12=rand(nx,nu);
% B21 = rand(np,nw); B22=rand(np,nu);
% D11 = rand(nz,nw); D12=rand(nz,nu);
% D21 = rand(ny,nw); D22=rand(ny,nu);
% Bx = rand(nrL1/2,nx); Bw = rand(nrL1/2,nw); Bu = rand(nrL1/2,nu);
% E6 = rand(np,nx); bC1 = rand(nz,nx); bC2 = rand(ny,nx);
% for i=1:N
%     ntmp=i*ni{i+1};
%     E1{i} = rand(nx,ntmp);E2{i} = rand(nx,ntmp);
%     E4{i} = rand(np,ntmp);E5{i} = rand(np,ntmp);
%     C1{i} = rand(nz,ntmp);C2{i} = rand(nz,ntmp);
%     C4{i} = rand(ny,ntmp);C5{i} = rand(ny,ntmp);
%     B1{i} = rand(nrL1/2,ntmp);B2{i} = rand(nrL1/2,ntmp);
% end
% 
% for i=0:N
%     ntmp=(i+1)*ni{i+1};
%     B3{i+1} = rand(nrL1/2,ntmp);
%     E3{i+1} = rand(nx,ntmp);
%     C3{i+1} = rand(nz,ntmp);
%     C6{i+1} = rand(ny,ntmp);
%     A0{i+1} = rand(np,ntmp);
%     A1{i+1} = rand(np,ntmp);
%     A2{i+1} = rand(np,ntmp);
% end

%% Heat equation - stable
% pvar s theta;
% ni = {0,0,1}; N = length(ni)-1;
% nz = 0; ny = 0; nw =0; nu=0; nx = 0;
% a = 0; b = 1;
% 
% np = 0; nrL1 = 0;
% for i =0:N
% np=np+ni{i+1};
% nrL1 = nrL1+i*ni{i+1};
% end
% nrL1 = 2*nrL1;
% 
% A = rand(nx);
% B11 = rand(nx,nw); B12=rand(nx,nu);
% B21 = rand(np,nw); B22=rand(np,nu);
% D11 = rand(nz,nw); D12=rand(nz,nu);
% D21 = rand(ny,nw); D22=rand(ny,nu);
% Bx = rand(nrL1/2,nx); Bw = rand(nrL1/2,nw); Bu = rand(nrL1/2,nu);
% E6 = rand(np,nx); bC1 = rand(nz,nx); bC2 = rand(ny,nx);
% for i=1:N
%     ntmp=i*ni{i+1};
%     E1{i} = zeros(nx,ntmp);E2{i} = zeros(nx,ntmp);
%     E4{i} = zeros(np,ntmp);E5{i} = zeros(np,ntmp);
%     C1{i} = zeros(nz,ntmp);C2{i} = zeros(nz,ntmp);
%     C4{i} = zeros(ny,ntmp);C5{i} = zeros(ny,ntmp);
%     B1{i} = zeros(nrL1/2,ntmp);B2{i} = zeros(nrL1/2,ntmp);
% end
% 
% for i=0:N
%     ntmp=(i+1)*ni{i+1};
%     B3{i+1} = zeros(nrL1/2,ntmp);
%     E3{i+1} = zeros(nx,ntmp);
%     C3{i+1} = zeros(nz,ntmp);
%     C6{i+1} = zeros(ny,ntmp);
%     A0{i+1} = zeros(np,ntmp);
%     A1{i+1} = zeros(np,ntmp);
%     A2{i+1} = zeros(np,ntmp);
% end
% A0{3} = [0 0 1];
% B1 = {[1;0],[0;0]}; B2 = {[0;1],[0;0]};
% 
% stability =1;
% stability_dual=1;
%% forced convection over solid mass - stable
% pvar s theta;
% ni = {0,0,1}; N = length(ni)-1;
% nz = 0; ny = 0; nw =0; nu=0; nx = 0;
% a = 0; b = 1;
% 
% np = 0; nrL1 = 0;
% for i =0:N
% np=np+ni{i+1};
% nrL1 = nrL1+i*ni{i+1};
% end
% nrL1 = 2*nrL1;
% 
% A = rand(nx);
% B11 = rand(nx,nw); B12=rand(nx,nu);
% B21 = rand(np,nw); B22=rand(np,nu);
% D11 = rand(nz,nw); D12=rand(nz,nu);
% D21 = rand(ny,nw); D22=rand(ny,nu);
% Bx = rand(nrL1/2,nx); Bw = rand(nrL1/2,nw); Bu = rand(nrL1/2,nu);
% E6 = rand(np,nx); bC1 = rand(nz,nx); bC2 = rand(ny,nx);
% for i=1:N
%     ntmp=i*ni{i+1};
%     E1{i} = zeros(nx,ntmp);E2{i} = rand(nx,ntmp);
%     E4{i} = zeros(np,ntmp);E5{i} = zeros(np,ntmp);
%     C1{i} = rand(nz,ntmp);C2{i} = rand(nz,ntmp);
%     C4{i} = rand(ny,ntmp);C5{i} = rand(ny,ntmp);
%     B1{i} = rand(nrL1/2,ntmp);B2{i} = rand(nrL1/2,ntmp);
% end
% 
% for i=0:N
%     ntmp=(i+1)*ni{i+1};
%     B3{i+1} = zeros(nrL1/2,ntmp);
%     E3{i+1} = rand(nx,ntmp);
%     C3{i+1} = rand(nz,ntmp);
%     C6{i+1} = rand(ny,ntmp);
%     A0{i+1} = rand(np,ntmp);
%     A1{i+1} = zeros(np,ntmp);
%     A2{i+1} = zeros(np,ntmp);
% end
% A0{3} = [0 0 1];
% B1 = {[1;0],[0;0]}; B2 = {[0;1],[0;1]};
% 
% stability =1;
%% Timoshenko beam Equation

%% Chemical Reactor - stable
% pvar s theta;
% ni = {0,0,1}; N = length(ni)-1;
% nz = 0; ny = 0; nw =0; nu=0; nx = 1;
% a = 0; b = 1;
% 
% np = 0; nrL1 = 0;
% for i =0:N
% np=np+ni{i+1};
% nrL1 = nrL1+i*ni{i+1};
% end
% nrL1 = 2*nrL1;
% 
% A = -0.2;
% B11 = rand(nx,nw); B12=rand(nx,nu);
% B21 = rand(np,nw); B22=rand(np,nu);
% D11 = rand(nz,nw); D12=rand(nz,nu);
% D21 = rand(ny,nw); D22=rand(ny,nu);
% Bx = zeros(nrL1/2,nx); Bw = zeros(nrL1/2,nw); Bu = zeros(nrL1/2,nu);
% E6 = zeros(np,nx); bC1 = zeros(nz,nx); bC2 = zeros(ny,nx);
% for i=1:N
%     ntmp=i*ni{i+1};
%     E1{i} = zeros(nx,ntmp);E2{i} = zeros(nx,ntmp);
%     E4{i} = zeros(np,ntmp);E5{i} = zeros(np,ntmp);
%     C1{i} = zeros(nz,ntmp);C2{i} = zeros(nz,ntmp);
%     C4{i} = zeros(ny,ntmp);C5{i} = zeros(ny,ntmp);
%     B1{i} = zeros(nrL1/2,ntmp);B2{i} = zeros(nrL1/2,ntmp);
% end
% 
% for i=0:N
%     ntmp=(i+1)*ni{i+1};
%     B3{i+1} = zeros(nrL1/2,ntmp);
%     E3{i+1} = zeros(nx,ntmp);
%     C3{i+1} = zeros(nz,ntmp);
%     C6{i+1} = zeros(ny,ntmp);
%     A0{i+1} = zeros(np,ntmp);
%     A1{i+1} = zeros(np,ntmp);
%     A2{i+1} = zeros(np,ntmp);
% end
% A0{1} = [-1]; A0{2}= [-1]; A0{3} = [0];
% B1 = {[1;0],[0;0]}; B2 = {[0;0],[0;1]}; E6 = 1;
% 
% stability =1;

%% Nonlocal boundary example 1 - stable
% pvar s theta;
% ni = {0,0,1}; N = length(ni)-1;
% nz = 0; ny = 0; nw =0; nu=0; nx = 0;
% a = 0; b = 1;
% 
% np = 0; nrL1 = 0;
% for i =0:N
% np=np+ni{i+1};
% nrL1 = nrL1+i*ni{i+1};
% end
% nrL1 = 2*nrL1;
% 
% A = [];
% B11 = rand(nx,nw); B12=rand(nx,nu);
% B21 = rand(np,nw); B22=rand(np,nu);
% D11 = rand(nz,nw); D12=rand(nz,nu);
% D21 = rand(ny,nw); D22=rand(ny,nu);
% Bx = zeros(nrL1/2,nx); Bw = zeros(nrL1/2,nw); Bu = zeros(nrL1/2,nu);
% E6 = zeros(np,nx); bC1 = zeros(nz,nx); bC2 = zeros(ny,nx);
% for i=1:N
%     ntmp=i*ni{i+1};
%     E1{i} = zeros(nx,ntmp);E2{i} = zeros(nx,ntmp);
%     E4{i} = zeros(np,ntmp);E5{i} = zeros(np,ntmp);
%     C1{i} = zeros(nz,ntmp);C2{i} = zeros(nz,ntmp);
%     C4{i} = zeros(ny,ntmp);C5{i} = zeros(ny,ntmp);
%     B1{i} = zeros(nrL1/2,ntmp);B2{i} = zeros(nrL1/2,ntmp);
% end
% 
% for i=0:N
%     ntmp=(i+1)*ni{i+1};
%     B3{i+1} = zeros(nrL1/2,ntmp);
%     E3{i+1} = zeros(nx,ntmp);
%     C3{i+1} = zeros(nz,ntmp);
%     C6{i+1} = zeros(ny,ntmp);
%     A0{i+1} = zeros(np,ntmp);
%     A1{i+1} = zeros(np,ntmp);
%     A2{i+1} = zeros(np,ntmp);
% end
% A0{1} = [0]; A0{2}= [0]; A0{3} = [1];
% B1 = {[1;0],[0;0]}; B2 = {[0;0],[0;1]}; E6 = [];
% B3 = {[0;-(s)^2],[0;0],[0;0]};
% 
% stability =1;

%% Nonlocal boundary example 2 - not stable?
pvar s theta;
ni = {0,0,1}; N = length(ni)-1;
nz = 0; ny = 0; nw =0; nu=0; nx = 0;
a = 0; b = 1;

np = 0; nrL1 = 0;
for i =0:N
np=np+ni{i+1};
nrL1 = nrL1+i*ni{i+1};
end
nrL1 = 2*nrL1;

A = [];
B11 = rand(nx,nw); B12=rand(nx,nu);
B21 = rand(np,nw); B22=rand(np,nu);
D11 = rand(nz,nw); D12=rand(nz,nu);
D21 = rand(ny,nw); D22=rand(ny,nu);
Bx = zeros(nrL1/2,nx); Bw = zeros(nrL1/2,nw); Bu = zeros(nrL1/2,nu);
E6 = zeros(np,nx); bC1 = zeros(nz,nx); bC2 = zeros(ny,nx);
for i=1:N
    ntmp=i*ni{i+1};
    E1{i} = zeros(nx,ntmp);E2{i} = zeros(nx,ntmp);
    E4{i} = zeros(np,ntmp);E5{i} = zeros(np,ntmp);
    C1{i} = zeros(nz,ntmp);C2{i} = zeros(nz,ntmp);
    C4{i} = zeros(ny,ntmp);C5{i} = zeros(ny,ntmp);
    B1{i} = zeros(nrL1/2,ntmp);B2{i} = zeros(nrL1/2,ntmp);
end

for i=0:N
    ntmp=(i+1)*ni{i+1};
    B3{i+1} = zeros(nrL1/2,ntmp);
    E3{i+1} = zeros(nx,ntmp);
    C3{i+1} = zeros(nz,ntmp);
    C6{i+1} = zeros(ny,ntmp);
    A0{i+1} = zeros(np,ntmp);
    A1{i+1} = zeros(np,ntmp);
    A2{i+1} = zeros(np,ntmp);
end
A0{1} = [0]; A0{2}= [0]; A0{3} = [1];
B1 = {[0;0],[0;0]}; B2 = {[0;0],[0;0]}; E6 = [];
B3 = {[1;(s)],[0;0],[0;0]};

stability =1;