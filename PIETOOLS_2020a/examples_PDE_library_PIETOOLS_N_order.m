clc; clear;
% Random ODE-PDE generator
pvar s theta;
ni = {1,2,3,5}; N = length(ni)-1;
nz = 1; ny = 2; nw =2; nu=1; nx = 2;
a = 0; b = 1;

np = 0; nrL1 = 0;
for i =0:N
np=np+ni{i+1};
nrL1 = nrL1+i*ni{i+1};
end
nrL1 = 2*nrL1;

A = rand(nx);
B11 = rand(nx,nw); B12=rand(nx,nu);
B21 = rand(np,nw); B22=rand(np,nu);
D11 = rand(nz,nw); D12=rand(nz,nu);
D21 = rand(ny,nw); D22=rand(ny,nu);
Bx = rand(nrL1/2,nx); Bw = rand(nrL1/2,nw); Bu = rand(nrL1/2,nu);
E6 = rand(np,nx); bC1 = rand(nz,nx); bC2 = rand(ny,nx);
for i=1:N
    ntmp=i*ni{i+1};
    E1{i} = rand(nx,ntmp);E2{i} = rand(nx,ntmp);
    E4{i} = rand(np,ntmp);E5{i} = rand(np,ntmp);
    C1{i} = rand(nz,ntmp);C2{i} = rand(nz,ntmp);
    C4{i} = rand(ny,ntmp);C5{i} = rand(ny,ntmp);
    B1{i} = rand(nrL1/2,ntmp);B2{i} = rand(nrL1/2,ntmp);
end

for i=0:N
    ntmp=(i+1)*ni{i+1};
    B3{i+1} = rand(nrL1/2,ntmp);
    E3{i+1} = rand(nx,ntmp);
    C3{i+1} = rand(nz,ntmp);
    C6{i+1} = rand(ny,ntmp);
    A0{i+1} = rand(np,ntmp);
    A1{i+1} = rand(np,ntmp);
    A2{i+1} = rand(np,ntmp);
end

convert_PIETOOLS_PDE_N_order;