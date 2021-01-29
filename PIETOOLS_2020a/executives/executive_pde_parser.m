function PDE_params = executive_pde_parser(PDE_text, n_ODE, n_PDE, n_dis, n_inp)

% remove whitespaces
PDE_text = regexprep(PDE_text,'\s+','');

% first identify the time derivative of the state on LHS
dynamics_state = regexp(PDE_text,'[a-zA-Z0-9]+t=','match');

% strip LHS to reduce regex complexity
PDE_text = regexprep(PDE_text,dynamics_state,'');

% strip 't='
dynamics_state = dynamics_state{1}(1:end-2);


if dynamics_state(1)=='x' % equation is an ODE
    n_number = str2double(dynamics_state(2:end)); %find state number to arrange the matrices
    if n_number>n_ODE
        error('State number exceeds the number of ODE states');
    end
    ode_coeffs = regexp(PDE_text,'[+-]?[\w*^()/]*x[0-9]+','match'); %ode terms
    pde_bcoeffs = regexp(PDE_text,'[+-]?[\w*^()/]*X[0-9]+(\(a\)|\(b\))','match'); %pde boundary terms
    dis_coeffs = regexp(PDE_text,'[+-]?[\w*^()/]*w[0-9]+','match'); %disturbance terms
    inp_coeffs = regexp(PDE_text,'[+-]?[\w*^()/]*u[0-9]+','match'); %input terms
    
    
    pvar s; 
    % extract coeffs
    for i=1:n_ODE
        token = regexp(ode_coeffs,['[+-]?[\w*^()]*x' num2str(i)], 'match'); 
        token = token(~cellfun('isempty',token)); 
        token(cellfun('isempty',token))=[];% assuming no repeats
        if length(token)==1
            cff = regexprep(token{1}{:},['*x' num2str(i)],'');
            A(i) = str2num(cff);
        else
            A(i)=0;
        end
    end
    
    for i=1:n_dis
        token = regexp(ode_coeffs,['[+-]?[\w*^()]*w' num2str(i)], 'match'); 
        token = token(~cellfun('isempty',token)); 
        token(cellfun('isempty',token))=[];% assuming no repeats
        if length(token)==1
            cff = regexprep(token{1}{:},['*w' num2str(i)],'');
            B11(i) = str2num(cff);
        else
            B11(i)=0;
        end
    end
    
    % extract coeffs
    for i=1:n_ODE
        token = regexp(ode_coeffs,['[+-]?[\w*^()]*u' num2str(i)], 'match'); 
        token = token(~cellfun('isempty',token)); 
        token(cellfun('isempty',token))=[];% assuming no repeats
        if length(token)==1
            cff = regexprep(token{1}{:},['*u' num2str(i)],'');
            B12(i) = str2num(cff);
        else
            B12(i)=0;
        end
    end
    
elseif dynamics_state(1)=='X' % equation is an PDE
    n_number = str2double(dynamics_state(2:end)); %find state number to arrange the matrices
    if n_number>n_PDE
        error('State number exceeds the number of PDE states');
    end
    ode_coeffs = regexp(PDE_text,'[+-]?[\w*^()/]*x[0-9]+','match'); %ode terms
    pde_bcoeffs = regexp(PDE_text,'[+-]?[\w*^()/]*X[0-9]+(\(a\)|\(b\))','match'); %pde boundary terms
    dis_coeffs = regexp(PDE_text,'[+-]?[\w*^()/]*w[0-9]+','match'); %disturbance terms
    inp_coeffs = regexp(PDE_text,'[+-]?[\w*^()/]*u[0-9]+','match'); %input terms
    pde_coeffs = regexp(PDE_text,'[+-]?[\w*^()/]*X[0-9]+\(s\)','match'); %pde boundary terms
else 
    error('Incorrect input format for the ODE-PDE');
end
end