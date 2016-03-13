%function [idname,nodes,id,cutpoint]=eval(Tree,X,subtrees)
function [idname nodes]=eval(Tree,X,subtrees)
%EVAL Compute fitted value for decision tree applied to data.
%   YFIT = EVAL(T,X) takes a classification or regression tree T and a
%   matrix X of predictor values, and produces a vector YFIT of predicted
%   response values. For a regression tree, YFIT(J) is the fitted response
%   value for a point having the predictor values X(J,:).  For a
%   classification tree, YFIT(J) is the class into which the tree would
%   assign the point with data X(J,:).
%
%   YFIT = EVAL(T,X,SUBTREES) takes an additional vector SUBTREES of
%   pruning levels, with 0 representing the full, unpruned tree.  T must
%   include a pruning sequence as created by the CLASSREGTREE constructor or
%   the PRUNE method. If SUBTREES has K elements and X has N rows, then the
%   output YFIT is an N-by-K matrix, with the Ith column containing the
%   fitted values produced by the SUBTREES(I) subtree.  SUBTREES must be
%   sorted in ascending order. (To compute fitted values for a tree that is
%   not part of the optimal pruning sequence, first use PRUNE to prune the
%   tree.)
%
%   [YFIT,NODE] = EVAL(...) also returns an array NODE of the same size
%   as YFIT containing the node number assigned to each row of X.  The
%   VIEW method can display the node numbers for any node you select.



[nr,nc] = size(X);

if nargin<3
   subtrees = 0;
elseif prod(size(subtrees))>length(subtrees)
   error('stats:treeval:BadSubtrees','SUBTREES must be a vector.');
elseif any(diff(subtrees)<0)
   error('stats:treeval:BadSubtrees','SUBTREES must be sorted.');
end

if ~isempty(Tree.prunelist)
   prunelist = Tree.prunelist;
elseif ~isequal(subtrees,0)
   Tree = prune(Tree);
else
   prunelist = repmat(Inf,size(Tree.node));
end

ntrees = length(subtrees);
toProcess.nodes       = zeros(nr,ntrees);

nnodes = length(Tree.node);
toProcess.rows        = cell(nnodes,1);
toProcess.thisnode    = zeros(nnodes,1);
toProcess.endcol      = zeros(nnodes,1);

toProcess.rows{1}     = 1:nr;
toProcess.thisnode(1) = 1;
toProcess.endcol(1)   = ntrees;
toProcess.iend        = 1;


iproc = 0;
while iproc < toProcess.iend
    iproc = iproc + 1;
    toProcess = doapply(Tree,X,toProcess,iproc,subtrees,prunelist);
end
nodes = toProcess.nodes;
nodes1=nodes(find(nodes~=0));
nodes2=nodes(find(nodes==0));
id(find(nodes~=0)) = Tree.class(nodes1);
id(find(nodes==0))=-1;


idname(find(nodes~=0)) = Tree.classname(id(find(nodes~=0)));

if ~isempty(nodes2)
    ind=find(nodes==0);
    for j=1:length(ind)
      idname(ind(j)) = {'-1'} ;
    end
end




%------------------------------------------------
function toProcess = doapply(Tree,X,toProcess,iproc,subtrees,prunelist)
%DOAPPLY Apply classification rule to specified rows starting at a node
%           given by iproc.
%   X, PRUNELIST, and SUBTREES do not change.
%
%   toProcess is a structure that specifies the process order and 
%       contains results:
%   NODES has one row per observation and one column per subtree.
%   ROWS describes the subset of X and NODES to consider. 
%   1:ENDCOL are columns of NODES and the elements of SUBTREES to consider.

rows     = toProcess.rows{iproc};
thisnode = toProcess.thisnode(iproc);
endcol   = toProcess.endcol(iproc);
iend     = toProcess.iend;

splitvar      = Tree.var(thisnode);
cutoff        = Tree.cut{thisnode};
kids          = Tree.children(thisnode,:);
prunelevel    = prunelist(thisnode);
categ_cols= Tree.catcols;
categ_grp=Tree.catgrp;

if ~isempty(categ_cols)
    categ_len=length(categ_cols)/categ_grp;
end
% For how many of the remaining trees is this a terminal node?
if splitvar==0      % all, if it's terminal on the unpruned tree

    ncols = endcol;
else                % some, if it's terminal only after pruning
   ncols = sum(subtrees(1:endcol) >= prunelevel);
end
if ncols>0          % for those trees, assign the node level now
   toProcess.nodes(rows,(endcol-ncols+1:endcol)) = thisnode;
  
   endcol = endcol - ncols;
end

% Now deal with non-terminal nodes
if endcol > 0
   % Determine if this point goes left, goes right, or stays here
     
   if splitvar>0     % continuous variable
       x = X(rows,abs(splitvar));
      isleft = (x < cutoff);
      isright = ~isleft;
      ismissing = isnan(x);
      on_line=[];
   end
   if  splitvar<0     % categorical variable
     cutpt=cutoff{1};
     
      value=[];
      cols=(-splitvar-1)*categ_len+categ_cols(1):-splitvar*categ_len+categ_cols(1)-1;
      
      
    
       x=[X(rows,cols) ones(length(rows),1)];
       
       value=cutpt*x';
    
      
      isleft = find(value<0);
      isright = find(value>0);
      
      
   end


 
 
   % Left children
   
   subrows = rows(isleft);
   if ~isempty(subrows)
       iend = iend + 1;
       toProcess.rows{iend}     = subrows;
       toProcess.thisnode(iend) = kids(1);
       toProcess.endcol(iend)   = endcol;
   end
   
   % Right children
   subrows = rows(isright );
   if ~isempty(subrows)
       iend = iend + 1;
       toProcess.rows{iend}     = subrows;
       toProcess.thisnode(iend) = kids(2);
       toProcess.endcol(iend)   = endcol;
   end
   
   
   % Update the last filled index
   toProcess.iend = iend;
end