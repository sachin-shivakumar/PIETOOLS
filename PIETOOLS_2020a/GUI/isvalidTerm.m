function logval = isvalidTerm(exp, termtype, eqntype, multiplierterm, intvar)
% This function tests whether a mathematical expression exp is a valid
% term to be added to dynamics or BCs
% if logval is 0, then exp is a valid polynomial
% logval is 2, exp has unidentified polynomial variable different from s or theta
% logval is 3, exp has polynomial variable with non-integer power
% logval is 4, exp has other error such as mismatched brackets, incorrect
% operands etc.

logval=0;
pvar s theta;

exp = convertCharsToStrings(exp); % convert to string

if strcmp(eqntype,'PDE')
    logval=0;
    if strcmp(termtype,'Multiplier')&&~contains(exp,'(\theta)')...
            &&~contains(multiplierterm,'(\theta)')...
            &&~contains(exp,'theta')
        tmp=0;
    elseif strcmp(termtype,'Integral')
        if (contains(exp,'theta')||contains(multiplierterm,'(\theta)')) && strcmp(intvar,'ds')
            logval=2;
        end
    else
        logval=2;
    end
    
elseif strcmp(eqntype, 'ODE') || strcmp(eqntype, 'BC')
    if strcmp(termtype,'Multiplier')&&~contains(multiplierterm,'(s)')...
            &&~contains(multiplierterm,'(\theta)')
        try
            tmp = double(eval(exp));
        catch err
            logval=2;
        end
    elseif strcmp(termtype,'Integral')
        if (contains(exp,'s')||contains(multiplierterm,'(s)'))&&strcmp(intvar,'d\theta')
            logval=2;
        end
        if (contains(exp,'theta')||contains(multiplierterm,'(\theta)'))&&strcmp(intvar,'ds')
            logval=2;
        end
    else 
        logval=2;
    end
end

end