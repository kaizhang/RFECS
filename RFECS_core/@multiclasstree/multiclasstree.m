classdef multiclasstree
%multiclasstree Create a classification and regression tree object.
%   T = multiclasstree(X,Y) creates a decision tree T for predicting response Y as
%   a function of predictors X.  X is an N-by-M matrix of predictor values.
%   If Y is a vector of N response values, then multiclasstree performs
%   regression.  If Y is a categorical variable, character array, or cell
%   array of strings, multiclasstree performs classification.  Either way, T is
%   a binary tree where each non-terminal node is split based on the values
%   of a column of X.  NaN values in X or Y are taken to be missing values,
%   and observations with any missing values are not used in the fit.
%
%   T = multiclasstree(X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
%
%   For all trees:
%      'categorical' Vector of indices of the columns of X that are to be
%                   treated as unordered categorical variables
%    
%      'names'      A cell array of names for the predictor variables,
%                   in the order in which they appear in the X matrix
%                   from which the tree was created
%      'prune'      'on' (default) to compute the full tree and the optimal
%                   sequence of pruned subtrees, or 'off' for the full tree
%                   without pruning
%      'minparent'  A number K such that impure nodes must have K or more
%                   observations to be split (default 10).
%      'minleaf'    A minimal number of observations per tree leaf
%                   (default=1). If both 'minparent' and 'minleaf' are
%                   supplied, the setting which results in larger leaves is
%                   used: MINPARENT = MAX(MINPARENT,2*MINLEAF) 
%      'nvartosample' Number of predictor variables randomly selected
%                     for each split. By default all variables are
%                     considered for each decision split.
%      'mergeleaves' 'on' (default) to merge leaves that originate from the
%                    same parent node and give the sum of risk values
%                    greater or equal to the risk associated with the
%                    parent node. If 'off', leaves are not merged.
%
%   For regression trees only:
%      'qetoler'     Defines tolerance on quadratic error per node for
%                    regression trees. Splitting nodes stops when quadratic
%                    error per node drops below QETOLER*QED, where QED is
%                    the quadratic error for the entire data computed
%                    before the decision tree is grown: QED = NORM(Y-YBAR)
%                    with YBAR estimated as the average of the input array
%                    Y. Default = 1e-6.  
%
%   For classification trees only:
%      'cost'       Square matrix C, where C(i,j) is the cost of classifying
%                   a point into class j if its true class is i (default
%                   has C(i,j)=1 if i~=j, and C(i,j)=0 if i=j).  Alternatively,
%                   this value can be a structure S having two fields:  S.group
%                   containing the group names as a categorical variable,
%                   character array, or cell array of strings; and S.cost
%                   containing the cost matrix C.
%      
%      'priorprob'  Prior probabilities for each class, specified as a
%                   vector (one value for each distinct group name) or as a
%                   structure S with two fields:  S.group containing the group
%                   names as a categorical variable, character array, or cell
%                   array of strings; and S.prob containing a vector of
%                   corresponding probabilities.
%

    properties(SetAccess='protected',GetAccess='protected')
      
          node = zeros(0,1);
        parent = zeros(0,1);
         class = zeros(0,1);
           var = zeros(0,1);
           cut = cell(0,1);
      children = zeros(0,2);
      nodeprob = zeros(0,1);
       nodeerr = zeros(0,1);
      nodesize = zeros(0,1);
         npred = 0;
       catcols = [];
       catgrp = [];
        prior = [];
      nclasses = 1;
          cost = [];
     classprob = [];
      impurity = [];
    classcount = [];
     classname = {};
     prunelist = zeros(0,1);
         alpha = [];
    ntermnodes = [];
         names = {};
prunecriterion = '';
rows           = [];
    end
    
    methods
        function a= multiclasstree(x,y,varargin)
            if nargin==1 && isa(x,'struct')
                a= struct2tree(a,x);            % convert from struct to tree
            else
                error(nargchk(2,Inf,nargin,'struct'));
                a = treefit_new(a,x,y,varargin{:});  % calls local version of treefit_new
            end
        end % multiclasstree constructor
        
    end % methods block

    methods(Hidden = true)
        function b = properties(a)
            b = {'node'; 'class'; ...
                 'npred'; 'catcols'; 'catgrp';'prior'; 'nclasses'; 'cost'; ...
                 'classname'; 'prunelist'; 'alpha'; 'ntermnodes'; ...
                 'names'; 'impurity'; 'prunecriterion';};
        end
        function b = fieldnames(a)
            b = properties(a);
        end
        
        % Methods that we inherit, but do not want
        function a = fields(varargin),     throwUndefinedError(); end
        function a = ctranspose(varargin),  throwUndefinedError(); end
        function a = transpose(varargin),  throwUndefinedError(); end
        function a = permute(varargin),   throwUndefinedError(); end
        function a = reshape(varargin),   throwUndefinedError(); end
        function a = cat(varargin),        throwNoCatError(); end
        function a = horzcat(varargin),    throwNoCatError(); end
        function a = vertcat(varargin),    throwNoCatError(); end
    end
    methods(Hidden = true, Static = true)
        function a = empty(varargin)
            error(['stats:' mfilename ':NoEmptyAllowed'], ...
                  'Creation of empty %s objects is not allowed.',upper(mfilename));
        end
    end

    methods(Static=true,Hidden=true)
        function [X,Y,classnames,Yorig] = preparedata2(X,Y)
            %preparedata2 Perform initial data checks and transformations.
            %
            % [XOUT,YOUT] = preparedata2(X,Y,DOCLASS) processes matrix of predictor
            % values X and vector of true response Y. For regression, Y is a numeric
            % vector of regressed values. For classification, Y is a categorical vector
            % of true class labels represented as an array of either char's or integers
            % or a cell array of strings. X is a numeric matrix with one row per
            % observation and one column per input variable. The number of rows in X
            % must match the number of elements in Y. DOCLASS is a logical flag, true
            % for classification and false for regression.
            %
            % preparedata2 removes rows with nothing but missing values (NaN's) in X,
            % rows without a valid class label in Y (for classification), and rows with
            % NaN's in Y (for regression). For regression, it simply returns X and Y
            % with these rows removed as XOUT and YOUT.
            %
            % For classification, [XOUT,YOUT] = preparedata2(X,Y,true)
            % returns a numeric vector YOUT of group indices for groups
            % found in input labels Y.
            %
            % [XOUT,YOUT,CLASSNAMES,YORIG] = preparedata2(X,Y,true) also returns class
            % names and YORIG, categorical Y data for classification. The class names
            % are for the groups found in input labels Y. YORIG is always a cell arry
            % of strings with the original Y labels while YOUT is always a vector of
            % numeric indices irrespective of the type of Y.
            %
            % See also GRP2IDX.
        
            
            % Make class labels
            classnames = {};
            Yorig = Y;
         
                [Y,classnames] = grp2idx(Y);
                t = isnan(Y);
                if any(t)
                    Yorig(t) = [];
                    Y(t) = [];
                    X(t,:) = [];
                end
                if isempty(X)
                    error('stats:multiclasstree:preparedata2:NoData',...
                        'No data remaining in X and Y after removing NaN groups.');
                end
                if ischar(Yorig)
                    Yorig = cellstr(Yorig);
                end
          
        end
        
        function varnames = preparevars(varnames,nvars)
            %PREPAREVARS Perform initial variable check and transformations.
            %
            % VARNAMES = PREPAREVARS(VARNAMES,NVARS) checks the size of the input
            % vector of variable names VARNAMES against the expected size NVARS,
            % converts VARNAMES to a cell array if VARNAMES are char's and assigns
            % default variable names if none are supplied.
            if ~isempty(varnames)
                if ischar(varnames)
                    varnames = cellstr(varnames);
                end
                if ~iscellstr(varnames) || numel(varnames)~=nvars
                    error('stats:multiclasstree:preparevars:BadNames',...
                        'Variable names must be a character array or cell array with %d strings.',...
                        nvars);
                end
            else
                varnames = strcat('x',strread(sprintf('%d\n',1:nvars),'%s\n'));
            end
        end
    end
end % classdef

function throwNoCatError()
error(['stats:' mfilename ':NoCatAllowed'], ...
      'Concatenation of %s objects is not allowed.  Use a cell array to contain multiple objects.',upper(mfilename));
end

function throwUndefinedError()
st = dbstack;
name = regexp(st(2).name,'\.','split');
error(['stats:' mfilename ':UndefinedFunction'], ...
      'Undefined function or method ''%s'' for input arguments of type ''%s''.',name{2},mfilename);
end



function Tree=treefit_new(Tree,X,y,varargin)

% Process inputs
% Categorical is the components of the vector variable and
% categorical_group is the number of categorical variables

okargs =   {'priorprob'   'cost'  'splitcriterion' ...
            'splitmin' 'minparent' 'minleaf' ...
            'nvartosample' 'mergeleaves' 'categorical' 'categorical_group'  'prune'  ...
            'qetoler' 'names'};
defaults = {[]        []      'gdi'                        ...
            []        []        1                          ...
            'all'     'on'     []      []     'on'            ...
            1e-6      {}};
%[eid,emsg,Prior,Cost,Criterion,splitmin,minparent,minleaf,...
 %   nvartosample,Merge,categ,categ_grp,Prune,qetoler,names,extra] = ...
  %  dfswitchyard('statgetargs',okargs,defaults,varargin{:});

  %if ~isempty(emsg)
   % error(sprintf('stats:treefit_new:%s',eid),emsg);
%end
[Prior,Cost,Criterion,splitmin,minparent,minleaf,...
    nvartosample,Merge,categ,categ_grp,Prune,qetoler,names,extras] = ...
    internal.stats.parseArgs(okargs,defaults,varargin{:});


% Preprocess data
[X,y,cnames] = multiclasstree.preparedata2(X,y);
[N,nvars2] = size(X);
if(isempty(categ))
    nvars=nvars2;
else
    nvars=nvars2-length(categ)+categ_grp;
end

% Process variable names
names = multiclasstree.preparevars(names,nvars);

% Fill out criterion, class labels and matrix for classification

   switch(Criterion)
    %                Criterion function   Is it an impurity measure?
    %                ------------------   --------------------------
    case 'gdi',      critfun = @gdi;      isimpurity = 1;
    
    otherwise,     error('stats:treefit_new:BadSplitCriterion',...
                         'Bad value for ''splitcriterion'' parameter.')
   end
   
   % Get binary matrix, C(i,j)==1 means point i is in class j
   nclasses = length(cnames);
   C = false(N,nclasses);
   C(sub2ind([N nclasses],(1:N)',y)) = 1;
   Nj = sum(C,1);


% Check and adjust node sizes.
if ~isempty(splitmin) && ~isempty(minparent)
    error('stats:treefit_new:BadNodeSize',...
        'One cannot use splitmin and minparent at the same time.');
end
if ~isempty(splitmin)
    if ~isnumeric(splitmin)
        error('stats:treefit_new:BadNodeSize',...
            '''splitmin'' argument must be numeric.');
    end
    minparent = splitmin;
end
if ~isnumeric(minparent)
    error('stats:treefit_new:BadNodeSize',...
        '''minparent'' argument must be numeric.');
end
if isempty(minparent)
    minparent = 10;
end
if ~isnumeric(minleaf)
    error('stats:treefit_new:BadNodeSize',...
        '''minleaf'' argument must be numeric.');
end
if minleaf<1
    error('stats:treefit_new:BadNodeSize',...
        'Minimal size for leaf node must be greater or equal to 1.');
end
if minparent<1
    error('stats:treefit_new:BadNodeSize',...
        'Minimal size for parent node must be greater or equal to 1.');
end
minparent = max(minparent,2*minleaf);

% Set the number of vars to be selected at random for splits
% Get number of features to sample
success = false;
if strcmpi(nvartosample,'all')
    nvartosample = nvars;
    success = true;
end
if isnumeric(nvartosample) && nvartosample>0 && nvartosample<=nvars
    nvartosample = ceil(nvartosample);
    success = true;
end
if ~success
    error('stats:treefit_new:BadNumberOfRandomFeatures',...
        'Number of features to be selected at random for decision splits must be ''all'' or a number between 1 and %i: %s',...
        nvars,num2str(nvartosample));
end
nusevars = nvartosample;

% Check prune paremeter
if ~strcmpi(Prune,'on') && ~strcmpi(Prune,'off')
    error('stats:treefit_new:BadParamName',...
        '''prune'' parameter must be either ''on'' or ''off''.');
end

% Check merge paremeter
if ~strcmpi(Merge,'on') && ~strcmpi(Merge,'off')
    error('stats:treefit_new:BadParamName',...
        '''mergeleaves'' parameter must be either ''on'' or ''off''.');
end

% Tree structure fields ([C] only for classification trees):
%  .method     method
%  .node       node number
%  .parent     parent node number
%  .class      class assignment for points in this node if treated as a leaf
%  .var        column j of X matrix to be split, or 0 for a leaf node,
%              or -j to treat column j as categorical
%  .cut        cutoff value for split (Xj<cutoff goes to left child node),
%              or a cell array of left and right values if var is negative
%  .children   matrix of child nodes (2 cols, 1st is left child)
%  .nodeprob   probability p(t) for this node
%  .nodeerr    resubstitution error estimate r(t) for this node
%  .nodesize   number of points at this node
%  .prunelist  list of indices that define pruned subtrees.  One entry per
%              node.  If prunelist(j)=k then, at the kth level of pruning,
%              the jth node becomes a leaf (or drops off the tree if its
%              parent also gets pruned).
%  .alpha      vector of complexity parameters for each pruning cut
%  .ntermnodes vector of terminal node counts for each pruning cut
%  .classprob  [C] vector of class probabilities
%  .classname  [C] names of each class
%  .classcount [C] count of members of each class
%  .nclasses   [C] number of classes
%  .cost       [C] misclassification cost

M = 2*ceil(N/minleaf)-1;% number of tree nodes for space reservation
nodenumber = zeros(M,1);
parent = zeros(M,1);
yfitnode = zeros(M,1);
cutvar = zeros(M,1);
cutpoint = cell(M,1);
children = zeros(M,2);
nodeprob = zeros(M,1);
resuberr = zeros(M,1);
nodesize = zeros(M,1);

   classprob = zeros(M,nclasses);
   classcount = zeros(M,nclasses);
   if isimpurity==1
       impurity = zeros(M,1);
   end

iscat = zeros(nvars2,1); iscat(categ) = 1;

nodenumber(1) = 1;

assignednode = cell(M,1);% list of instances assigned to this node
assignednode{1} = 1:N;
nextunusednode = 2;


   Prior = Prior(:)';
   haveprior = true;
   if isempty(Prior)
      Prior = Nj / N;
      haveprior = false;
   elseif isequal(Prior,'equal')
      Prior = ones(1,nclasses) / nclasses;

   elseif isstruct(Prior)
      if ~isfield(Prior,'group') || ~isfield(Prior,'prob')
         error('stats:treefit_new:BadPrior',...
              'Missing field in structure value for ''priorprob'' parameter.');
      end
      idx = getclassindex(cnames,Prior.group);
      if any(idx==0)
         j = find(idx==0);
         error('stats:treefit_new:BadPrior',...
               'Missing prior probability for group ''%s''.',cnames{j(1)});
      end
      Prior = Prior.prob(idx);
   end
   if length(Prior)~=nclasses || any(Prior<0) || sum(Prior)==0 ...
                              || ~isnumeric(Prior)
      error('stats:treefit_new:BadPrior',...
            'Value of ''priorprob'' parameter must be a vector of %d probabilities.',...
            nclasses);
   elseif all(Prior==0 | sum(C,1)==0)
      error('stats:treefit_new:BadPrior',...
            'The ''priorprob'' parameter assigns all probability to unobserved classes.');
   else
      Prior = Prior / sum(Prior);
   end

   % Get default or specified misclassification costs
   havecosts = true;
   if isempty(Cost)
      Cost = ones(nclasses) - eye(nclasses);
      havecosts = false;
   else
      if isstruct(Cost)
         if ~isfield(Cost,'group') || ~isfield(Cost,'cost')
            error('stats:treefit_new:BadCost',...
                  'Missing field in structure value for ''cost'' parameter.');
         end
         idx = getclassindex(cnames,Cost.group);
         if any(idx==0)
            j = find(idx==0);
            error('stats:treefit_new:BadCost',...
                  'Missing misclassification cost for group ''%s''.',...
                          cnames{j(1)});
         end
         Cost = Cost.cost(idx,idx);
      end
      if ~isequal(size(Cost),nclasses*ones(1,2))
         error('stats:treefit_new:BadCost',...
               'Misclassification cost matrix must be %d-by-%d.',...
                       nclasses,nclasses);
      elseif any(diag(Cost)~=0)
         error('stats:treefit_new:BadCost',...
            'Misclassification cost matrix must have zeros on the diagonal.');
      elseif any(Cost<0)
         error('stats:treefit_new:BadCost',...
            'Misclassification cost matrix must contain non-negative values.');
      end
   end
   
   % Adjust priors if required to take misclassification costs into account
   adjprior = Prior;
   if havecosts
      Cj = sum(Cost,2)';
      pc = Cj .* Prior;
      adjprior = pc / sum(pc);
   end
   pratio = adjprior ./ max(1,Nj);


% Keep processing nodes until done
tnode = 1;
while(tnode < nextunusednode)
   % Record information about this node
   noderows = assignednode{tnode};
  
   Nnode = length(noderows);
   Cnode = C(noderows,:);
   
      % Compute class probabilities and related statistics for this node
      Njt = sum(Cnode,1);    % number in class j at node t
      Pjandt = Prior .* Njt ./ max(1,Nj);   %max() to get 0/1 if Njt=Nj=0
      Pjgivent = Pjandt / sum(Pjandt);
      misclasscost = Pjgivent * Cost;
      [mincost,nodeclass] = min(misclasscost);
      yfitnode(tnode) = nodeclass;
      Pt = sum(Pjandt);
      nodeprob(tnode) = Pt;
      classprob(tnode,:) = Pjgivent;
      classcount(tnode,:) = Njt;
      impure = sum(Pjgivent>0)>1;
      if isimpurity==1
          impurity(tnode) = feval(critfun,classprob(tnode,:));
      end
   
   bestcrit          = -Inf;
   nodesize(tnode)   = Nnode;
   resuberr(tnode)   = mincost;
   cutvar(tnode)     = 0;
   cutpoint{tnode}   = 0;
   children(tnode,:) = 0;
   
   % Consider splitting this node
   if (Nnode>=minparent) && impure      % split only large impure nodes
      Xnode = X(noderows,:);
      bestvar = 0;
      bestcut = 0;

      % Reduce the number of predictor vars as specified by nvarstosample
      varmap = 1:nvars;
      if nusevars < nvars
         varmap = randsample(nvars,nusevars);
      end
      
      % Find the best of all possible splits
      %Select all non-vector variables
   
      nusevars1 = intersect(find(iscat==0),varmap);
      nusevars2 = find(iscat==1);
      
      
      for ivar=1:length(nusevars1)
          % Index of variable to split on
          %jvar = varmap(ivar);
          
         jvar=nusevars1(ivar);
       
              
              idxnan = isnan(Xnode(:,jvar));
              idxnotnan = find(~idxnan);
              if isempty(idxnotnan)
                  continue;
              end
              [x,idxsort] = sort(Xnode(idxnotnan,jvar));
              idx = idxnotnan(idxsort);
              c = Cnode(idx,:);
              
              Ccum = cumsum(c,1);
              
              % Determine if there's anything to split along this variable
              maxeps = max(eps(x(1)), eps(x(end)));
              if x(1)+maxeps > x(end)
                  continue;
              end
              
              % Accept only splits on rows with distinct values
              rows = find( x(1:end-1) + ...
                  max([eps(x(1:end-1)) eps(x(2:end))],[],2) < x(2:end) );
              if isempty(rows)
                  continue;
              end
              
              
              
              % Reduce the list of splits to use only rows with enough class
              % counts. For categorical vars need to consider all split
              % permutations => reducing Ccum based on row counts doesn't work
              
              Ctot = sum(Ccum(end,:));% total counts from all classes
              Crow = sum(Ccum,2);% cumsum of class counts on each row
              rows = rows( Crow(rows)>=minleaf & (Ctot-Crow(rows))>=minleaf );
              if isempty(rows)
                  continue;
              end
         
         
         % For classification, keep only rows which have more than one
         % class between them
        
             nrow = length(rows);
             keep = false(nrow,1);
             % Check between 1st instance and 1st kept row
             if any(any(c(1:rows(1)-1,:) ~= c(2:rows(1),:),2))
                 keep(1) = true;
             end
             % Check all kept rows
             for ir=1:nrow-1
                 if any(any(c(rows(ir):rows(ir+1)-1,:) ~= ...
                         c(rows(ir)+1:rows(ir+1),:),2))
                     keep(ir:ir+1) = true;
                 end
             end
             % Check between last kept row and last instance
             Nx = length(x);
             if any(any(c(rows(end):Nx-1,:) ~= c(rows(end)+1:Nx,:),2))
                 keep(end) = true;
             end
             % Reduce the list of rows
             rows = rows(keep);
         
         
         % If there are missing values, recompute quantities
         %   needed for splits ignoring these
       
         ybar0 = [];
         
         Pt0 = Pt;
         
         if any(idxnan)
           
                 Njt0 = sum(c,1);
                 Pjandt0 = Prior .* Njt0 ./ max(1,Nj);
                 Pt0 = sum(Pjandt0);
            
         end
         
         % Compute numerical criteria for splits.
         % Use one routine for regression and classification.
         [critval,cutval]=RCcritval(x,c,Ccum,rows,pratio,Pt0,ybar0,isimpurity,critfun,bestcrit,minleaf);

         % Change best split if this one is best so far
         if critval>bestcrit
            bestcrit = critval;
            bestvar = jvar;
            bestcut = cutval;
         end
      end
      
        if(~isempty(categ_grp))
      categ_grp2=setdiff(varmap,find(iscat==0))-length(find(iscat==0));
      
      categ_len=length(nusevars2)/categ_grp;
      
      for j=1:length(categ_grp2)
          
          %check if the node can be split
         
          vector1=Xnode(:,nusevars2(categ_len*(categ_grp2(j)-1)+1:categ_len*categ_grp2(j)));
         
          
         Cnode2=zeros(size(Cnode,1),1);
          Cnode2(find(Cnode(:,1)==1))=1;
         Cnode2(find(Cnode(:,2)==1))=0;
          try
          
          [class_lda err postprob logp coeff ]=classify(vector1,vector1,Cnode2,'diaglinear');
          catch ME1
              idSegLast=regexp(ME1.identifier, '(?<=:)\w+$','match');
              if (strcmp(idSegLast,'stats:classify:BadVariance')==1)
                  continue;
                  
              end
          end
          
          cutval={[coeff(1,2).linear;coeff(1,2).const]'};
         
          
          %split on left node and right node
        Csplit1= [length(intersect(find(Cnode2==1),find((coeff(1,2).linear)'*vector1'+coeff(1,2).const <0))) length(intersect(find(Cnode2==0),find((coeff(1,2).linear)'*vector1'+coeff(1,2).const <0)))];
         Csplit2= [length(intersect(find(Cnode2==1),find((coeff(1,2).linear)'*vector1'+coeff(1,2).const >0))) length(intersect(find(Cnode2==0),find((coeff(1,2).linear)'*vector1'+coeff(1,2).const >0)))];
          
         
        P1 = pratio .* Csplit1;
        P2= pratio.*Csplit2;
        
        Ptleft  = sum(P1,2);
        Ptright = sum(P2,2);
       
        nclasses = size(P1,2);
        wuns = ones(1,nclasses);
        P1 = P1 ./ Ptleft(:,wuns);   %repmat(Ptleft,1,nclasses);
        P2 = P2./Ptright(:,wuns);
        % Get left/right node probabilities
        Pleft = Ptleft ./ Pt;
        Pright = Ptright./Pt;
        
        % Evaluate criterion as impurity or otherwise
        if isimpurity
            crit = - Pleft.*feval(critfun,P1)-Pright.*feval(critfun,P2);
           
           
             critval=crit;
        end
  
      if(critval>bestcrit)
        
         bestcrit = critval;
         bestvar2=categ_grp2(j);
         bestvar=nusevars2(categ_len*(categ_grp2(j)-1)+1:categ_len*categ_grp2(j));
         bestcut = cutval;
       leftside=find((coeff(1,2).linear)'*vector1'+coeff(1,2).const <0);
        rightside=find((coeff(1,2).linear)'*vector1'+coeff(1,2).const >0);
         
         
      end
      end
      end
 
          
          

      % Split this node using the best rule found
      if bestvar~=0
          
         x = Xnode(:,bestvar);
   
         if (length(find(iscat(bestvar)==1))==0)
            
            cutvar(tnode) = bestvar;
            leftside = x<=bestcut;
            rightside = ~leftside;
            
         else
            cutvar(tnode) = -bestvar2;          % negative indicates cat. var. split
            %leftside = ismember(x,bestcut{1});
            
            %rightside = ismember(x,bestcut{2});
         end
         cutpoint{tnode} = bestcut;
         children(tnode,:) = nextunusednode + (0:1);
         
         
         
         assignednode{nextunusednode} = noderows(leftside);
         assignednode{nextunusednode+1} = noderows(rightside);
         nodenumber(nextunusednode+(0:1)) = nextunusednode+(0:1)';
         parent(nextunusednode+(0:1)) = tnode;
         nextunusednode = nextunusednode+2;
      end
   end
   
   tnode = tnode + 1;
end

topnode        = nextunusednode - 1;

Tree.node      = nodenumber(1:topnode);
Tree.parent    = parent(1:topnode);
Tree.class     = yfitnode(1:topnode);
Tree.var       = cutvar(1:topnode);
Tree.cut       = cutpoint(1:topnode);
Tree.children  = children(1:topnode,:);
Tree.nodeprob  = nodeprob(1:topnode);
Tree.nodeerr   = resuberr(1:topnode);
Tree.nodesize  = nodesize(1:topnode);
Tree.npred     = nvars;
Tree.catcols   = categ;
Tree.catgrp = categ_grp;
Tree.names     = names;
Tree.rows = assignednode(1:topnode);
%Tree.rows  = assignednode(1:topnode);
   Tree.prior     = Prior;
   Tree.nclasses  = nclasses;
   Tree.cost      = Cost;
   Tree.classprob = classprob(1:topnode,:);
   Tree.classcount= classcount(1:topnode,:);
   Tree.classname = cnames;
   if isimpurity==1
       Tree.impurity = impurity(1:topnode);
   end


if isequal(Merge,'on')
    Tree = mergeleaves(Tree); % merge leaves with same class
end

if isequal(Prune,'on')        % compute optimal pruning sequence if requested
   Tree = prune(Tree);
end

end

%----------------------------------------------------
function v=gdi(p)
%GDI Gini diversity index

v=1-sum(p.^2,2);
end


%----------------------------------------------------
    function [critval,cutval] = ...
            RCcritval(x,c,Ccum,rows,pratio,Pt,ybar, ...
            isimpurity,critfun,bestcrit,minleaf)
        
        % How many splits?
        nsplits = length(rows);
        % Split between each pair of distinct ordered values
        Csplit1 = Ccum(rows,:);
        
        
        
        % Complementary splits
        Csplit2 = Ccum(size(Ccum,1)*ones(nsplits,1),:) - Csplit1;
        
      
        
        
        % Get left/right class probabilities or mean values at each split
        
        % Classification
        temp = pratio(ones(nsplits,1),:); %repmat(pratio,nsplits,1);
        P1 = temp .* Csplit1;
        P2 = temp .* Csplit2;
        Ptleft  = sum(P1,2);
        Ptright = sum(P2,2);
        nclasses = size(P1,2);
        wuns = ones(1,nclasses);
        P1 = P1 ./ Ptleft(:,wuns);   %repmat(Ptleft,1,nclasses);
        P2 = P2 ./ Ptright(:,wuns);  %repmat(Ptright,1,nclasses);
        
        % Get left/right node probabilities
        Pleft = Ptleft ./ Pt;
        Pright = 1 - Pleft;
        
        % Evaluate criterion as impurity or otherwise
        if isimpurity
            crit = - Pleft.*feval(critfun,P1);
            t = (crit>bestcrit); % compute 2nd term only if it would make a difference
            if any(t)
                crit(t) = crit(t) - Pright(t).*feval(critfun,P2(t,:));
            end
        else
            crit = feval(critfun, Pleft, P1, Pright, P2);
        end
        
        
        % Get best split point
        [critval,maxloc] = max(crit);
        if critval<bestcrit
            cutval = NaN;
            return;
        end
        
        % Get the cut value
        
        cutloc = rows(maxloc);
        cutval = (x(cutloc) + x(cutloc+1))/2;
        
    end

% --------------------------------------
    
function Tree = mergeleaves(Tree)
        %MERGELEAVES Merge leaves that originate from the same parent node and give
        % the sum of risk values greater or equal to the risk associated with the
        % parent node.
        
        N = length(Tree.node);
        isleaf = (Tree.var==0)';   % no split variable implies leaf node
        isntpruned = true(1,N);
        doprune = false(1,N);
        Risk = risk(Tree)';
        adjfactor = (1 - 100*eps(class(Risk)));
        
        % Work up from the bottom of the tree
        while(true)
            branches = find(~isleaf & isntpruned);
            twig = branches(sum(isleaf(Tree.children(branches,:)),2) == 2);
            if isempty(twig)
                break;            % must have just the root node left
            end
            
            % Find twigs to ''unsplit'' if the error of the twig is no larger
            % than the sum of the errors of the children
            Rtwig = Risk(twig);
            kids = Tree.children(twig,:);
            Rsplit = sum(Risk(kids),2);
            unsplit = Rsplit >= Rtwig'*adjfactor;
            if any(unsplit)
                % Mark children as pruned, and mark twig as now a leaf
                isntpruned(kids(unsplit,:)) = 0;
                twig = twig(unsplit);   % only these to be marked on next 2 lines
                isleaf(twig) = 1;
                doprune(twig) = 1;
            else
                break;
            end
        end
        
        % Remove splits that are useless
        if any(doprune)
            Tree = prune(Tree,'nodes',find(doprune));
        end
    end

% ------------------------------------
function idx = getclassindex(cnames,g)
%GETCLASSINDEX Find indices for class names in another list of names
%   IDX = GETCLASSINDEX(CNAMES,G) takes a list CNAMES of class names
%   (such as the grouping variable values in the treefit_new or classify
%   function) and another list G of group names (as might be supplied
%   in the ''prior'' argument to those functions), and finds the indices
%   of the CNAMES names in the G list.  CNAMES should be a cell array
%   of strings.  G can be numbers, a string array, or a cell array of
%   strings

% Convert to common string form, whether input is char, cell, or numeric
if isnumeric(g)
   g = cellstr(strjust(num2str(g(:)), 'left'));
elseif ~iscell(g)
   g = cellstr(g);
end

nclasses = length(cnames);
idx = zeros(1,nclasses);

% Look up each class in the grouping variable.
for i = 1:nclasses
   j = strmatch(cnames{i}, g, 'exact');
   if ~isempty(j)
      idx(i) = j(1);
   end
end
end

% ------------------------------------
function treeobj=struct2tree(treeobj,S)
% Copy fields from structure S to tree object treeobj

% Look at all fields required for regression or classification trees
allfields = {   'node'     'parent'   'class'   'var' ...
             'cut'      'children' 'nodeprob' 'nodeerr' ...
             'nodesize' 'npred'    'catcols' 'catgrp'  ...
             'nclasses' 'prior'    'cost'     ...
             'classprob' 'classcount' 'classname'};
fn = fieldnames(S);

    nrequired = numel(allfields);

for j=1:nrequired
    fname = allfields{j};
    treeobj.(fname) = S.(fname);
end

% Look at optional fields
optionalfields = {'names' 'prunelist' 'alpha' 'ntermnodes' ...
    'impurity' 'prunecriterion'};
for j=1:numel(optionalfields)
    fname = optionalfields{j};
    if isfield(S,fname)
        treeobj.(fname) = S.(fname);
    end
end
end





