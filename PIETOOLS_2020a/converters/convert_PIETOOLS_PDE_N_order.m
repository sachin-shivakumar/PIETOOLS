%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_PDE_MMP.m     PIETOOLS 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_PDE_MMP is an alternative version of the PDE converter 
% file which performs the following two tasks.
% 1) It verifies the dimension compatibility of input parameters of ODE-PDE
% and sets any missing parameters to zero.
% 2) It converts the input ODE-PDE representation to a PIE
% representation. 
%
% A Partial Integral Equation is defined by 11 PI operators as
%
% BT1op \dot{w}(t)+BT2op \dot{u}(t)+Top \dot{x}(t)= Aop x(t) + B1op u(t)+ + B2op w(t)
%                                             z(t)= C1op x(t) + D11op u(t)+ + D12op w(t)
%                                             y(t)= C2op x(t) + D21op u(t)+ + D22op w(t)
%
% This script takes a user-defined PDE system in the format outlined in the
% header of solver_PIETOOLS_PDE and converts it to a PIE by defining the 11
% PI operators {BT1op,BT2op,Top,Aop,B1op,B2op,C1op,D11op,D12op,C2op,D21op,D22op} 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following script performs error checking operations and sets
% undefined operators to empty opbjects.
% initialize_PIETOOLS_PDE


% Converts ODE-PDE to PIE and defines PI operators
disp('Converting ODE-PDE to PIE');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define auxiliary variables to convert the ODE-PDE to a PIE
X=[a b];

zwu=zeros(nz,nu);zuw=zwu';
zwx=zeros(nz,nx);zxw=zwx';
Iw=eye(nw);

np = 0; nrL1 = 0; nrL2 =0;
for i =0:N
np=np+ni{i+1};
nrL1 = nrL1+2*i*ni{i+1};
nrL2 = nrL2+(i+1)*ni{i+1};
end

ncL1=np;
ncL2=np;
 
%%%%%%%%%%%%%%%%%%%% Defining Primal Dynamics %%%%%%%%%%%%%%%%%%%%%%%
% Rewriting in PI notation for ODE-PDE
bE1_a=[]; bE1_b =[]; bE3_a=[];  bE3_b=[]; 
bC10_a=[]; bC10_b=[]; bC20_a=[]; bC20_b=[]; 
bC11=[];bC12=[]; bE2=[]; bA0=[]; bA1=[]; bA2=[];
for i=0:N-1
    bE1_a =[bE1_a,E1{i+1}];bE1_b =[bE1_b,E2{i+1}];
    bE3_a =[bE3_a,E4{i+1}];bE3_b =[bE3_b,E5{i+1}];
    bC10_a =[bC10_a,C1{i+1}];bC10_b =[bC10_b,C2{i+1}];
    bC20_a =[bC20_a,C4{i+1}];bC20_b =[bC20_b,C5{i+1}];
end
bE1 = [bE1_a bE1_b]; bE3 = [bE3_a bE3_b];
bC10 = [bC10_a bC10_b]; bC20 = [bC20_a bC20_b];
for i=0:N
    bE2 =[bE2, E3{i+1}];
    bC11=[bC11, C3{i+1}];
    bC12=[bC12, C6{i+1}];
    bA0=[bA0, A0{i+1}];
    bA1=[bA1, A1{i+1}];
    bA2=[bA2, A2{i+1}];
end

Pb=[D11 D12 bC1 bC10;
    D21 D22 bC2 bC20;
    B11 B12 A  bE1];
Q1b=[bC11;
     bC12;
     bE2];
Q2b=[B21 B22 E6 bE3];
R0b=bA0;
R1b=bA1; R2b=bA2;

%%% At this point, we can construct the primal dynamics opvar
opvar Apop;
Apop.dim = [nz+ny+nx,nw+nu+nx+nrL1;np,nrL2]; Apop.var1 = s; Apop.var2 = theta; Apop.I = X;
Apop.P=Pb;
Apop.Q1=Q1b;
Apop.Q2=Q2b;
Apop.R.R0=R0b;
Apop.R.R1=R1b;
Apop.R.R2=R2b;


%%%%%%%%%%%%%%%%%%%% Converting PhiN X  %%%%%%%%%%%%%%%%%%%%%%%
bB1_a=[]; bB1_b=[]; bB2=[]; 

for i=0:N-1
    bB1_a =[bB1_a,B1{i+1}];bB1_b =[bB1_b,B2{i+1}];
end
bB1 = [bB1_a bB1_b];
for i=0:N
    bB2 =[bB2,B3{i+1}];
end


%%% The following auxiliary matrices are used in this section 
for k=1:N
    bT{k} =[];
    bQ{k} =[];
    for i=k-1:-1:0
        bTtemp=[];
        for j=0:i
            bTtemp = [bTtemp (s-a)^(j)/factorial(j)*eye(ni{k+1})];
        end
    bTtemp = [zeros(ni{k+1},(k-i-1)*ni{k+1}) bTtemp];
    bT{k} = [bT{k}; bTtemp];
    bQ{k} = [bQ{k}; (s-theta)^(i)/factorial(i)*eye(ni{k+1})];
    end
end

T21 =[]; Q21 =[];
for k=1:N
T21 = blkdiag(T21,bT{k});
Q21 = blkdiag(Q21,bQ{k});
end
Q21 = [zeros(nrL1/2,ni{1}),Q21];


T11 = [subs(T21,s,a);subs(T21,s,b)];
Q11 = [subs(Q21,s,a);subs(Q21,s,b)];


U2 = zeros(nrL2,np); U1 = eye(nrL2);

for i=N:-1:0
    indy=0;
    for j=0:i
        indy = indy+(j+1)*ni{j+1};
    end
    U1(:,indy-ni{i+1}+1:indy) =[];
end
for j=0:N
    indx=0;indy=0;
    for i=0:j
        indx = indx+(i+1)*ni{i+1};
        indy = indy+ni{i+1};
    end
    U2(indx-ni{j+1}+1:indx,indy-ni{j+1}+1:indy) = eye(ni{j+1});
end

B_T = double(bB1*T11+int(bB2*U1*T21,s,a,b));
pvar eta;
Q_T = bB1*Q11+subs(int(subs(bB2*U1*Q21,s,eta),eta,s,b),theta,s)+bB2*U2;

K=zeros(ni{1});
L0=eye(ni{1});
L1=zeros(ni{1});


for i=1:N
    bT1 = bT{i}; bQ1 = bQ{i};
    K = blkdiag(K,bT1(1:ni{i+1},:));
    L1 = blkdiag(L1,bQ1(1:ni{i+1},:));
    L0 = blkdiag(L0,zeros(ni{i+1}));
end
K = K(:,ni{1}+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for well-posedness of Boundary conditions

if rcond(B_T)<1e-15
    error('Defined boundary conditions are rank deficient or have prohibited boundary conditions. Your system is likely ill-posed.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
BTinv=inv(B_T);
Btemp=BTinv*[Bw Bu Bx];

J0 = K*BTinv*Bx;
J1 = K*BTinv*Bw;
J2 = K*BTinv*Bu;
G0 = L0; G2 = -K*BTinv*bB1*Q11; G1 = G2+L1;

J3 = U1*T21*BTinv*Bx;
J4 = U1*T21*BTinv*Bw;
J5 = U1*T21*BTinv*Bu;
H0 = U2; H2 = -U1*T21*BTinv*subs(Q_T,s,theta); H1 = H2+U1*Q21;

%%%%%%%%%%%%%%%% Construct the LHS opvar 
opvar Top TBw TBu;
Top.dim = [nx nx; np np]; Top.var1 = s; Top.var2 = theta; Top.I = X;
Top.P = eye(nx);
Top.Q2 = J0; Top.R.R0 = G0; Top.R.R1 = G1; Top.R.R2 = G2;

TBw.dim = [nx nw; np 0]; TBw.var1 = s; TBw.var2 = theta; TBw.I = X;
TBw.Q2 = J1; 

TBu.dim = [nx nu; np 0]; TBu.var1 = s; TBu.var2 = theta; TBu.I = X;
TBu.Q2 = J2; 


%%%%%%%%%%% Now construct the RHS opvar

opvar TDop TDBw TDBu;
TDop.dim = [nx nx; np np]; TDop.var1 = s; TDop.var2 = theta; TDop.I = X;
TDop.P = eye(nx);
TDop.Q2 = J3; TDop.R.R0 = H0; TDop.R.R1 = H1; TDop.R.R2 = H2;

TDBw.dim = [nx nw; np 0]; TDBw.var1 = s; TDBw.var2 = theta; TDBw.I = X;
TDBw.Q2 = J4; 

TDBu.dim = [nx nu; np 0]; TDBu.var1 = s; TDBu.var2 = theta; TDBu.I = X;
TDBu.Q2 = J5; 

Btemp=BTinv*[Bw Bu Bx];

 
Phf=[eye(nw+nu+nx);
    T11*Btemp];
Q1hf=[zeros(nw+nu+nx,np); (subs(subs(Q11,s,0),theta,s)-T11*BTinv*Q_T)] ;
Q2hf=[J4 J5 J3];
R0hf=H0;
R1hf=H1;
R2hf=H2;

%%% We now construct the Phfop
opvar Phfop;
Phfop.dim = [nw+nu+nx+nrL1, nw+nu+nx;nrL2, np]; Phfop.var1 = s; Phfop.var2 = theta; Phfop.I = X;
Phfop.P=Phf;
Phfop.Q1=Q1hf;
Phfop.Q2=Q2hf;
Phfop.R.R0=R0hf;
Phfop.R.R1=R1hf;
Phfop.R.R2=R2hf;

%%% We now construct the Ptop
Ptop=Apop*Phfop;



%%% Now we have to partition Ptop to get the desired pieces
D11op=Ptop(1:nz,1:nw);
D21op=Ptop((nz+1):(nz+ny),1:nw);
B1op=Ptop((nz+ny+1):(nz+ny+nx+np),1:nw);
D12op=Ptop(1:nz,(nw+1):(nw+nu));
D22op=Ptop((nz+1):(nz+ny),(nw+1):(nw+nu));
B2op=Ptop((nz+ny+1):(nz+ny+nx+np),(nw+1):(nw+nu));
C1op=Ptop(1:nz,(nw+nu+1):(nw+nu+nx+np));
C2op=Ptop((nz+1):(nz+ny),(nw+nu+1):(nw+nu+nx+np));
Aop=Ptop((nz+ny+1):(nz+ny+nx+np),(nw+nu+1):(nw+nu+nx+np));
TB1op=TBw;
TB2op=TBu;



nx1 = nx; nx2 = np;



%remove temporary opvars
clear Apop Phfop;