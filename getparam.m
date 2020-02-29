function [params, paramtable] = getparam(m1)
% uses model object, and returns a table of its reactions and parameters.

% because there can be multiple parameters of the same name, but in
% different reactions, we need to be careful and label the parameters
% by the reaction. Because the TXTL modeling toolbox does not name its
% reactions (it should!), we will use the entire reaction as the label.


% First, make a list of reaction scoped parameters 
count = 0;
for i = 1:length(m1.reactions)
    
    rx = get(m1.reactions(i), 'reaction'); %reaction string
    
    p = m1.reactions(i).KineticLaw.getparameters ;
    for j = 1:length(p)
        count = count+1;
        params{count,1} = get(p(j), 'Name');
        params{count,2} = get(p(j), 'Value');
        params{count,3} = rx;
        params{count,4} = m1.reactions(i);
        params{count,5} = p(j);
        
        rObj(count) = m1.reactions(i);
        pNum(count) = count;
        rName{count} = rx;
        pObj(count) = p(j);
        pName{count} = get(p(j), 'Name');
        pVal(count) = get(p(j), 'Value');
    end
end

% then add global parameters which are not already in some reaction
% to this list (April 4: why would global params (like 'TXTL_clpx_deg' and
% 'AGTPdeg_F' be in reaction scope too? Im confused. 
for i = 1:length(m1.parameters)
    idx = find(strcmp(m1.parameters(i).Name, params(:,1)),1);
    
    if isempty(idx)
        count = count+1;
        params{count,1} = m1.parameters(i).Name;
        params{count,2} = m1.parameters(i).value;
        params{count,3} = 'global';
        params{count,4} = '';
        params{count,5} = m1.parameters(i);
        
        rObj(count) = NaN; % create a reaction null -> null (based on the type of rest of the column. wow!)
        pNum(count) = count;
        rName{count} = 'global';
        pObj(count) = m1.parameters(i);
        pName{count} = m1.parameters(i).Name;
        pVal(count) = m1.parameters(i).value;
       
    end
end

% finally, compile a table array, to be output witht he standard cell
% array. Tables are often convenient. 
[pNum, pName, pVal, rName, rObj, pObj ] = deal(pNum', pName', pVal', rName', rObj', pObj' );
paramtable = table(pNum, pName, pVal, rName, rObj, pObj );




end
