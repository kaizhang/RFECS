function disp(t)
%DISP Display a CLASSREGTREE object.
%   DISP(T) prints the CLASSREGTREE object T.
%
%   See also CLASSREGTREE, CLASSREGTREE/VIEW.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.6.3 $  $Date: 2008/12/01 08:10:22 $

isLoose = strcmp(get(0,'FormatSpacing'),'loose');
maxWidth = get(0,'CommandWindowSize'); maxWidth = maxWidth(1);

% disp(struct(t));

if (isLoose), fprintf('\n'); end

% Get some information about the whole tree
maxnode = numel(t.node);
nd = 1 + floor(log10(maxnode)); % number of digits for node number
names = t.names;
if isempty(names)
    names = strcat('x',strread(sprintf('%d\n',1:t.npred),'%s\n'));
end



% Display information about each node
for j=1:maxnode
    if any(t.children(j,:))
        % branch node
        vnum = t.var(j);
        vname = names{abs(vnum)};
        cut = t.cut{j};
        kids = t.children(j,:);
        word=strcat(num2str(cut{1}(1)),',');
        for k=1:length(cut{1})
                word=strcat(word,',',num2str(cut{1}(k)));
        end
            cond = sprintf('%s<%s',vname,word);
            fprintf('%*d  if %s then node %d else node %d\n',nd,j,cond,kids);
        
    else
        % terminal node, display fit (regression) or class assignment
        
            fprintf('%*d  class = %s\n',nd,j,t.classname{t.class(j)});
  
    end
end
if (isLoose), fprintf('\n'); end