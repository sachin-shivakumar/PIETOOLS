%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examples_DDF_library_PIETOOLS.m     PIETOOLS 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This library file contains the definition of some common DDF systems
% drawn from the literature. To use one of the examples, simply uncomment
% the variable definitions for that example. Please make sure no other
% example is uncommented, or you will have problems.
%
% The examples are grouped into various problem types. These types are as
% follows:
% 1. Stability Analysis Tests
% 2. Hinf-gain Analysis
% 3. Optimal Control Problems
% 4. Optimal Estimator Design Problems
%
% The examples are typically called using a line in the main PIETOOLS_DDF.m
% file. Simply save your changes to the library file, uncomment the line
% examples_DDF_library_PIETOOLS.m in PIETOOLS_DDF.m and run PIETOOLS_DDF.m.
% Of course, you can also run the library file directly.
%
% Most examples in each problem type include the option to solve the LPI
% according to their problem type. Hence, it is recommended that you also
% uncomment the toggle as well. e.g. stability=1, Hinf_gain=1, etc.
%
% When relevant, we also include citations for each example, indicating the
% sources. The bibtex for each citation is included at the end of the
% library file.
%
% If you wish to include a new example in our library, please send it to us
% and we will include it in the next release. Please also include the
% citation information, if available.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STABILITY TEST EXAMPLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
<<<<<<< Updated upstream

% A pure Difference Equation r_1(t)=Drv{1}*Cv{1}*r_1(t-tau(1))+Drv{2}*Cv{2}*r_2(t-tau(2))
stability=1;
Drv{1}=.5; Cv{1}=1;Drv{2}=.25; Cv{2}=1;
tau(1)=1; tau(2)=2;
=======
% 
% % A pure Difference Equation r_1(t)=Drv{1}*Cv{1}*r_1(t-tau(1))+Drv{2}*Cv{2}*r_2(t-tau(2))
% stability=1;
% Drv{1}=.5; Cv{1}=1;Drv{2}=.25; Cv{2}=1;
% tau(1)=1; tau(2)=2;
% % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% An interface for declaring neutral-type systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A NDS Example Problem 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A0=[-2 .2 -.3 0 -.4;
%     .2 -3.8 0 .7 0;
%     .8 0 -1.6 0 0;
%     0 .8 -.6 -2 .3;
%     -1 -.1 -1.5 0 -1.8];
% Ai{1}=[-2.2 0 0 1 0;
%     1.6 -2.2 1.6 0 0;
%     -0.2 -0.2 -0.2 -0.2 -0.2;
%     0 0.4 -1.4 -3.4 1;
%     -0.2 0.4 -0.1 -1.1 -3.3];
% Ei{1}=[ 0.40888 0.00888 0.20888 -0.09112 -0.29112;
%     0 0.2 0 0 0.6;
%     -0.1 -0.4 0 -0.8 0;
%     0 0 -0.1 0 0;
%     0 0 0 -0.2 -0.1];
% 
% A0=[-.9 .2;.1 -.9];
% Ai{1}=[-1.1 -.2;-.1 -1.1];
% Ei{1}=[-.2 0;.2 -.1];
% tau(1)=1.29;
% 
% %%% stable for tau<1.3 - light settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Another NDS Example Problem 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A0=[-4 -1;0 -3];
% Ai{1}=[2 0;1 1];
% Ei{1}=+[.2 0;-.1 -.2];
% tau(1)=.5;
% 
% % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Another NDS Example Problem 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ki=10;kp=10;
% % a=.4;b=50;h=.2;d=.8;sig=.3;
% % al1=d+kp;gam1=b*ki*d+a*ki;
% % al2=(d-kp)*exp(sig*h);gam2=(b*ki*d-a*ki)*exp(sig*h);
% % beta1=(b*kp+a)*d+b*d^2+a*kp+ki;
% % beta2=((b*kp+a)*d-b*d^2-a*kp-ki)*exp(sig*h);
% % 
% % A0=(1/al1)*[0 al1;-sig^2*al1+sig*beta1-gam1 -beta1+2*sig*al1];
% % Ai{1}=(1/al1)*[0 0;-sig^2*al2+sig*beta2-gam2 -beta2+2*sig*al2];
% % Ei{1}=[0 0;0 -al2/al1];
% % tau(1)=h;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% A Neutral-Type System: 
% %%% \dot x(t)- \sum_i En{i}\dot x(t-\tau(i))=An x(t)+\sum_i Ani{i}x(t-\tau_i)
% %A0=[1];En{1}=[.6];Ani{1}=.2;
% ndim=size(An,1);
% stability=1;
% Cr{1}=[eye(ndim);An];
% A0=An;
% Bv=[Ani{1} En{1}];
% Cv{1}=eye(2*ndim);
% Drvi{1}=[zeros(ndim,2*ndim);Bv];
% % 

% %%% A 2-delay Neutral-Type System: 
% %%% \dot x(t)- \sum_i Ei{i}\dot x(t-\tau(i))=A0 x(t)+\sum_i Ai{i}x(t-\tau_i)
% %A0=[1];Ei{1}=[.6];Ai{1}=.2;
% Ei{1}=.2;Ei{2}=0;
% A0=-1;Ai{1}=0;Ai{2}=-2;
% h=.4; % max at .4448?
% tau=[h 2*h];

A0=[-2 0;0 -.9];Ai{1}=zeros(2);Ai{2}=[-1 0;-1 -1];
Ei{1}=[.1 0;0 .1];Ei{2}=zeros(2);
h=1; %max at 1.5804?
tau=[h 3*h];

%ndim=size(An,1);
stability=1;

initialize_PIETOOLS_NDS;
convert_PIETOOLS_NDS2DDF;

%minimize_PIETOOLS_DDF;

%convert_PIETOOLS_DDF;

% Cr{1}=[eye(ndim);An];
% Cr{2}=[eye(ndim);An];
% A=An;
% Bv=[Ani{1} En{1} Ani{2} En{2}];
% Cv{1}=[eye(2*ndim);zeros(2*ndim)];
% Cv{2}=[zeros(2*ndim);eye(2*ndim)];
% Drvi{1}=[zeros(ndim,4*ndim);Bv];
% Drvi{2}=[zeros(ndim,4*ndim);Bv];
>>>>>>> Stashed changes
% 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROLLER SYNTHESIS EXAMPLES
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %%% Multiple Showering People - tracking with integral control
% % % This is the DDF implementation
% Hinf_control=1
% n=5; %This is the number of users and can be changed
% shower_ex=2;
% nndelay=n;
% % the delay for user i is tau(i)=i;
% for i=1:nndelay
%     tau(i)=n;
% end
% % set all alphas to 1
% alpha=1;
% alphav=alpha*ones(n,1);
% % set all gammas to 1/n
% gamma=1/n;
% gammaM=gamma*ones(n,n);
% zn=zeros(n);
% 
% A0=[zn eye(n);zn zn];
% 
% Gam=zn;
% for i=1:n
%     for j=1:n
%         if i~=j
%             Gam(i,j)=alphav(j)*gammaM(i,j);
%         else
%             Gam(i,j)=-alphav(i);
%         end
%         
%     end
% end
% B1=[-eye(n);-Gam];
% B2=[zn;eye(n)];
% 
% %There are 3 option for the number of outputs in this example:
% %%%% Full outputs, Full disturbances
% % B1=[-eye(n);-Gam];
% % C1=[eye(n) zn;zn zn];
% % D11=[zn;zn];
% % D12=[zn;.1*eye(n)];
% 
% %%%%% 2 outputs, Full disturbances
% B1=[-eye(n);-Gam];
% C1=[ones(1,n) zeros(1,n);zeros(1,2*n)];
% D11=[zeros(2,n)];
% D12=[zeros(1,n);.1*ones(1,n)];
% % N=5 - light - gam = .714, IPM=17.4
% % N=30 - extreme - gam = 5.37, IPM=35,620
% 
% % Now for new DDF terms
% Bv=[zn;Gam];
% %D1v=zn
% 
% for i=1:nndelay
%     Cr{i}=[zeros(1,n+i-1) 1 zeros(1,n-i)];
% % Br1i=Br2i=Drvi=0;
% Cv{i}=zeros(n,1);Cv{i}(i,1)=1;
% %Cvdi=0;
% end
%%%gam_guess; % min in [.32 .38]
