function PDE_params = executive_pde_parser(PDE_Text, BC_Text, no, der_pde, nw, nu, nz, ny)

% remove whitespaces
PDE_Text = regexprep(PDE_Text,'\s+','');
BC_Text = regexprep(BC_Text,'\s+','');

der_pde = cellfun(@str2double, der_pde);

%initialize PDE parameters
init_PDE_params;

for k=1:size(PDE_Text,1)
% first identify the time derivative of the state on LHS
dynamics_state = regexp(PDE_Text{k},'[a-zA-Z0-9_]+=','match');
% strip LHS to reduce regex complexity
PDE_text = regexprep(PDE_Text{k},dynamics_state,'');
% strip '_t='
if strcmp(dynamics_state{1}(end-2:end),'_t=')
    dynamics_state = dynamics_state{1}(1:end-3);
else
    dynamics_state = dynamics_state{1}(1:end-1);
end

n_number = str2double(dynamics_state(2:end)); %find state number to arrange the matrices
ode_coeffs = regexp(PDE_text,'[+-]?[\w*^()+-]?*x[0-9]+','match'); %ode terms
pde_bcoeffs = regexp(PDE_text,'[+-]?[\w*^()/+-]?*X[0-9]+(_{s*})?(\(0\)|\(1\))','match'); %pde boundary terms
pde_icoeffs = regexp(PDE_text,'[+-]?(\\int_)+[\w*^()/+-]*X[0-9]+(_{s*})?([(s)]|[(Theta)])+([ds]|[dTheta])+','match'); %pde integral terms
dis_coeffs = regexp(PDE_text,'[+-]?[\w*^()/+-]?*w[0-9]+','match'); %disturbance terms
inp_coeffs = regexp(PDE_text,'[+-]?[\w*^()/+-]?*u[0-9]+','match'); %input terms
pvar s theta; 

if dynamics_state(1)=='x' % equation is an ODE
    % extract coeffs
    for i=1:length(ode_coeffs)
       col_num = str2double(regexprep(ode_coeffs{i},'[+-]?[\w*^()]*x',''));
       PDE.dyn.A(n_number,col_num) = str2double(regexprep(ode_coeffs{i},['*x' num2str(col_num)],''));
    end
    
    for i=1:length(dis_coeffs)
       col_num = str2double(regexprep(dis_coeffs{i},'[+-]?[\w*^()]*w',''));
       PDE.dyn.Bxw(n_number,col_num) = str2double(regexprep(dis_coeffs{i},['*w' num2str(col_num)],''));
    end
    
    for i=1:length(inp_coeffs)
       col_num = str2double(regexprep(inp_coeffs{i},'[+-]?[\w*^()]*u',''));
       PDE.dyn.Bxu(n_number,col_num) = str2double(regexprep(inp_coeffs{i},['*u' num2str(col_num)],''));
    end
    
    for i=1:length(pde_bcoeffs)
       state_num = str2double(regexprep(pde_bcoeffs{i},'([+-]?.**X)|(_{s*})|(\(0\))|(\(1\))',''));
       location_num = pde_bcoeffs{i}(end-2:end);
       derivative_num = regexp(pde_bcoeffs{i},'(_{s*})?','match');
       col_num = getcolnum(state_num, location_num, derivative_num, der_pde);
       PDE.dyn.Bxb(n_number,col_num) = str2double(regexprep(pde_bcoeffs{i},'(*X[0-9]+)|(_{s*})|(\(0\))|(\(1\))',''));
    end
    
    for i=1:length(pde_icoeffs)
       state_num = str2double(regexprep(pde_icoeffs{i},...
           '([+-]?.**X)|(_{s*})|(\(s\))|(\(Theta\))|(ds)|(dTheta)',''));
       location_num = 0;
       derivative_num = regexp(pde_icoeffs{i},'(_{s*})?','match');
       col_num = getcolnum2(state_num, location_num, derivative_num, der_pde);
       if isempty(derivative_num)
          derivative_num=0;
       else
          derivative_num = length(derivative_num{1})-3;
       end
       tmp = PDE.dyn.Bxr{derivative_num+1};
       tmp(n_number,col_num) = eval(strrep(regexprep(pde_icoeffs{i},...
           '(\\int_.\^.)|(*X[0-9](_{s*})?)|(\(s\)|\(Theta\))|((ds)|(dTheta))',''),'Theta','s'));
       PDE.dyn.Bxr{derivative_num+1} = tmp;
    end
    
elseif dynamics_state(1)=='X' % equation is an PDE
    Bvx = PDE.dyn.Dv(sum(PDE.n.np),no);
    Bvw = PDE.dyn.Dv(sum(PDE.n.np),nw);
    Bvu = PDE.dyn.Dv(sum(PDE.n.np),nu);
    % extract coeffs
    for i=1:length(ode_coeffs)
       col_num = str2double(regexprep(ode_coeffs{i},'[+-]?[\w*^()]*x',''));
       Bvx(n_number,col_num) = eval(regexprep(ode_coeffs{i},['*x' num2str(col_num)],''));
    end
    
    for i=1:length(dis_coeffs)
       col_num = str2double(regexprep(dis_coeffs{i},'[+-]?[\w*^()]*w',''));
       Bvw(n_number,col_num) = eval(regexprep(dis_coeffs{i},['*w' num2str(col_num)],''));
    end
    
    for i=1:length(inp_coeffs)
       col_num = str2double(regexprep(inp_coeffs{i},'[+-]?[\w*^()]*u',''));
       Bvu(n_number,col_num) = eval(regexprep(inp_coeffs{i},['*u' num2str(col_num)],''));
    end
    PDE.dyn.Dv = [Bvx Bvw Bvu];
    
    for i=1:length(pde_bcoeffs)
       state_num = str2double(regexprep(pde_bcoeffs{i},'([+-]?.**X)|(_{s*})|(\(0\))|(\(1\))',''));
       location_num = pde_bcoeffs{i}(end-2:end);
       derivative_num = regexp(pde_bcoeffs{i},'(_{s*})?','match');
       col_num = getcolnum(state_num, location_num, derivative_num, der_pde);
       PDE.dyn.Db(n_number,col_num) = eval(regexprep(pde_bcoeffs{i},'(*X[0-9]+)|(_{s*})|(\(0\))|(\(1\))',''));
    end
    
    pde_coeffs = regexp(PDE_text,'[+-]?[\w*^()/]*X[0-9]+(_s*)?\(s\)','match'); %pde terms
    
    for i=1:length(pde_coeffs)
       state_num = str2double(regexprep(pde_coeffs{i},...
           '([+-]?.**X)|(_{s*})|(\(s\))|(\(Theta\))',''));
       location_num = 0;
       derivative_num = regexp(pde_coeffs{i},'(_{s*})','match');
       col_num = getcolnum2(state_num, location_num, derivative_num, der_pde);
       if isempty(derivative_num)
          derivative_num=0;
       else
          derivative_num = length(derivative_num{1})-3;
       end
       tmp = PDE.dyn.A0{derivative_num+1};
       tmp(n_number,col_num) = eval(regexprep(pde_coeffs{i},['*X' num2str(col_num) '(_s*)?\(s\)'],''));
       PDE.dyn.A0{derivative_num+1}=tmp;
    end

    for i=1:length(pde_icoeffs)
       state_num = str2double(regexprep(pde_icoeffs{i},...
           '([+-]?.**X)|(_{s*})|(\(s\))|(\(Theta\))|(ds)|(dTheta)',''));
       location_num = 0;
       derivative_num = regexp(pde_icoeffs{i},'(_{s*})?','match');
       col_num = getcolnum2(state_num, location_num, derivative_num, der_pde);
       if isempty(derivative_num)
          derivative_num=0;
       else
          derivative_num = length(derivative_num{1})-3;
       end
       tmp = PDE.dyn.A1{derivative_num+1};
       tmp2 = PDE.dyn.A2{derivative_num+1};
       
       
       if contains(pde_icoeffs{i},'\int_0^s')
           poly_exp = strrep(regexprep(pde_icoeffs{i},...
           '(\\int_0\^s)|(*X[0-9](_{s*})?)|(\(s\)|\(Theta\))|((ds)|(dTheta))',''),'Theta','theta');
           tmp(n_number,col_num) = eval(poly_exp);
       elseif contains(pde_icoeffs{i},'\int_s^1')
           poly_exp2 = strrep(regexprep(pde_icoeffs{i},...
           '(\\int_s\^1)|(*X[0-9](_{s*})?)|(\(s\)|\(Theta\))|((ds)|(dTheta))',''),'Theta','theta');
           tmp2(n_number,col_num) = eval(poly_exp2);
       else 
           poly_exp3 = strrep(regexprep(pde_icoeffs{i},...
           '(\\int_0\^1)|(*X[0-9](_{s*})?)|(\(s\)|\(Theta\))|((ds)|(dTheta))',''),'Theta','theta');
           tmp(n_number,col_num) = eval(poly_exp3);
           tmp2(n_number,col_num) = eval(poly_exp3);
       end
       PDE.dyn.A1{derivative_num+1}=tmp;
       PDE.dyn.A2{derivative_num+1}=tmp2;
    end
    
elseif dynamics_state(1)=='z' %equation is a regulated output
    % extract coeffs
    for i=1:length(ode_coeffs)
       col_num = str2double(regexprep(ode_coeffs{i},'[+-]?[\w*^()]*x',''));
       PDE.out.C1(n_number,col_num) = str2double(regexprep(ode_coeffs{i},['*x' num2str(col_num)],''));
    end
    
    for i=1:length(dis_coeffs)
       col_num = str2double(regexprep(dis_coeffs{i},'[+-]?[\w*^()]*w',''));
       PDE.out.D1w(n_number,col_num) = str2double(regexprep(dis_coeffs{i},['*w' num2str(col_num)],''));
    end
    
    for i=1:length(inp_coeffs)
       col_num = str2double(regexprep(inp_coeffs{i},'[+-]?[\w*^()]*u',''));
       PDE.out.D1u(n_number,col_num) = str2double(regexprep(inp_coeffs{i},['*u' num2str(col_num)],''));
    end
    
    for i=1:length(pde_bcoeffs)
       state_num = str2double(regexprep(pde_bcoeffs{i},'([+-]?.**X)|(_{s*})|(\(0\))|(\(1\))',''));
       location_num = pde_bcoeffs{i}(end-2:end);
       derivative_num = regexp(pde_bcoeffs{i},'(_{s*})?','match');
       col_num = getcolnum(state_num, location_num, derivative_num, der_pde);
       PDE.out.Dzb(n_number,col_num) = str2double(regexprep(pde_bcoeffs{i},'(*X[0-9]+)|(_{s*})|(\(0\))|(\(1\))',''));
    end
    
    for i=1:length(pde_icoeffs)
       state_num = str2double(regexprep(pde_icoeffs{i},...
           '([+-]?.**X)|(_{s*})|(\(s\))|(\(Theta\))|(ds)|(dTheta)',''));
       location_num = 0;
       derivative_num = regexp(pde_icoeffs{i},'(_{s*})?','match');
       col_num = getcolnum2(state_num, location_num, derivative_num, der_pde);
       if isempty(derivative_num)
          derivative_num=0;
       else
          derivative_num = length(derivative_num{1})-3;
       end
       tmp = PDE.dyn.Dzr{derivative_num+1};
       tmp(n_number,col_num) = eval(strrep(regexprep(pde_icoeffs{i},...
           '(\\int_.\^.)|(*X[0-9](_{s*})?)|(\(s\)|\(Theta\))|((ds)|(dTheta))',''),'Theta','s'));
       PDE.out.Dzr{derivative_num+1}=tmp;
    end
else %equation is an observed output
    % extract coeffs
    for i=1:length(ode_coeffs)
       col_num = str2double(regexprep(ode_coeffs{i},'[+-]?[\w*^()]*x',''));
       PDE.out.C2(n_number,col_num) = str2double(regexprep(ode_coeffs{i},['*x' num2str(col_num)],''));
    end
    
    for i=1:length(dis_coeffs)
       col_num = str2double(regexprep(dis_coeffs{i},'[+-]?[\w*^()]*w',''));
       PDE.out.D2w(n_number,col_num) = str2double(regexprep(dis_coeffs{i},['*w' num2str(col_num)],''));
    end
    
    for i=1:length(inp_coeffs)
       col_num = str2double(regexprep(inp_coeffs{i},'[+-]?[\w*^()]*u',''));
       PDE.out.D2u(n_number,col_num) = str2double(regexprep(inp_coeffs{i},['*u' num2str(col_num)],''));
    end
    
    for i=1:length(pde_bcoeffs)
       state_num = str2double(regexprep(pde_bcoeffs{i},'([+-]?.**X)|(_{s*})|(\(0\))|(\(1\))',''));
       location_num = pde_bcoeffs{i}(end-2:end);
       derivative_num = regexp(pde_bcoeffs{i},'(_{s*})?','match');
       col_num = getcolnum(state_num, location_num, derivative_num, der_pde);
       PDE.out.Dyb(n_number,col_num) = str2double(regexprep(pde_bcoeffs{i},'(*X[0-9]+)|(_{s*})|(\(0\))|(\(1\))',''));
    end
    
    for i=1:length(pde_icoeffs)
       state_num = str2double(regexprep(pde_icoeffs{i},...
           '([+-]?.**X)|(_{s*})|(\(s\))|(\(Theta\))|(ds)|(dTheta)',''));
       location_num = 0;
       derivative_num = regexp(pde_icoeffs{i},'(_{s*})?','match');
       col_num = getcolnum2(state_num, location_num, derivative_num, der_pde);
       if isempty(derivative_num)
          derivative_num=0;
       else
          derivative_num = length(derivative_num{1})-3;
       end
       tmp = PDE.dyn.Dyr{derivative_num+1};
       tmp(n_number,col_num) = eval(strrep(regexprep(pde_icoeffs{i},...
           '(\\int_.\^.)|(*X[0-9](_{s*})?)|(\(s\)|\(Theta\))|((ds)|(dTheta))',''),'Theta','s'));
       PDE.dyn.Dyr{derivative_num+1} = tmp;
    end
end
end

for k=1:size(BC_Text,1)
% strip LHS to reduce regex complexity
BC_text = BC_Text{k}
% strip '0='
BC_text = BC_text(3:end);

ode_coeffs = regexp(BC_text,'[+-]?[\w*^()/]*x[0-9]+','match'); %ode terms
pde_bcoeffs = regexp(BC_text,'[+-]?[\w*^()/]*X[0-9]+(_{s*})?(\(0\)|\(1\))','match'); %pde boundary terms
pde_icoeffs = regexp(BC_text,'[+-]?(\\int_)+[\w*^()/]*X[0-9]+(_{s*})?([(s)]|[(Theta)])+([ds]|[dTheta])+','match'); %pde integral terms
dis_coeffs = regexp(BC_text,'[+-]?[\w*^()/]*w[0-9]+','match'); %disturbance terms
inp_coeffs = regexp(BC_text,'[+-]?[\w*^()/]*u[0-9]+','match'); %input terms
pvar s theta;
for i=1:length(ode_coeffs)
   col_num = str2double(regexprep(ode_coeffs{i},'[+-]?[\w*^()]*x',''));
   PDE.BC.Bx(k,col_num) = str2double(regexprep(ode_coeffs{i},['*x' num2str(col_num)],''));
end

for i=1:length(dis_coeffs)
   col_num = str2double(regexprep(dis_coeffs{i},'[+-]?[\w*^()]*w',''));
   PDE.BC.Bw(k,col_num) = str2double(regexprep(dis_coeffs{i},['*w' num2str(col_num)],''));
end

for i=1:length(inp_coeffs)
   col_num = str2double(regexprep(inp_coeffs{i},'[+-]?[\w*^()]*u',''));
   PDE.BC.Bu(k,col_num) = str2double(regexprep(inp_coeffs{i},['*u' num2str(col_num)],''));
end

for i=1:length(pde_bcoeffs)
   state_num = str2double(regexprep(pde_bcoeffs{i},'([+-]?.**X)|(_{s*})|(\(0\))|(\(1\))',''));
   location_num = pde_bcoeffs{i}(end-2:end);
   derivative_num = regexp(pde_bcoeffs{i},'(_{s*})?','match');
   col_num = getcolnum(state_num, location_num, derivative_num, der_pde);
   PDE.BC.B(k,col_num) = str2double(regexprep(pde_bcoeffs{i},'(*X[0-9]+)|(_{s*})|(\(0\))|(\(1\))',''));
end

for i=1:length(pde_icoeffs)
   state_num = str2double(regexprep(pde_icoeffs{i},...
       '([+-]?.**X)|(_{s*})|(\(s\))|(\(Theta\))|(ds)|(dTheta)',''));
   location_num = 0;
   derivative_num = regexp(pde_icoeffs{i},'(_{s*})?','match');
   col_num = getcolnum2(state_num, location_num, derivative_num, der_pde);
   if isempty(derivative_num)
      derivative_num=0;
   else
      derivative_num = length(derivative_num{1})-3;
   end
   tmp = PDE.BC.Br{derivative_num+1};
   tmp(k,col_num) = eval(strrep(regexprep(pde_icoeffs{i},...
       '(\\int_.\^.)|(*X[0-9](_{s*})?)|(\(s\)|\(Theta\))|((ds)|(dTheta))',''),'Theta','s'));
   PDE.BC.Br{derivative_num+1}=tmp;
end

end

PDE_params = PDE;
end
function col_loc = getcolnum(state_num, loc_num, deri_num, der_pde)
    if isempty(deri_num)
        deri_num=0;
    else
        deri_num = length(deri_num{1})-3;
    end
    
    loc_num = str2double(loc_num(2));
    
    for i=0:max(der_pde)
        n_pde(i+1) = sum(der_pde==i);
    end
    totalcols = sum(der_pde);
    der_pde_subset = der_pde(1:state_num); 
    der_pde_subset = der_pde_subset(der_pde_subset>=deri_num & der_pde_subset<=der_pde(state_num));
    col_loc = sum([1:deri_num-1].*n_pde(2:deri_num))+length(der_pde_subset);
    col_loc = col_loc+loc_num*totalcols;
end
function col_loc = getcolnum2(state_num, loc_num, deri_num, der_pde)
    if isempty(deri_num)
        deri_num=0;
    else
        deri_num = length(deri_num)-3;
    end
    for i=0:max(der_pde)
        n_pde(i+1) = sum(der_pde==i);
    end
    totalcols = sum(der_pde);
    der_pde_subset = der_pde(1:state_num); 
    der_pde_subset = der_pde_subset(der_pde_subset>=deri_num & der_pde_subset<=der_pde(state_num));
    col_loc = length(der_pde_subset);
end