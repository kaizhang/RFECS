function [cost,secost,ntnodes,bestlevel] = test_multi(Tree,TorCorR,X,Y,varargin)
%TEST Compute error rate for tree.
%   COST = TEST(T,'resubstitution') computes the cost of the tree T
%   using a resubstitution method.  T is a decision tree as created by
%   the CLASSREGTREE constructor.  The cost of the tree is the sum over all
%   terminal nodes of the estimated probability of that node times the
%   node's cost.  If T is a classification tree, the cost of a node is
%   the sum of the misclassification costs of the observations in
%   that node.  If T is a regression tree, the cost of a node is the
%   average squared error over the observations in that node.  COST is
%   a vector of cost values for each subtree in the optimal pruning
%   sequence for T.  The resubstitution cost is based on the same
%   sample that was used to create the original tree, so it under-
%   estimates the likely cost of applying the tree to new data.
%
%   COST = TEST(T,'test',X,Y) uses the predictor matrix X and
%   response Y as a test sample, applies the decision tree T to that
%   sample, and returns a vector COST of cost values computed for the
%   test sample.  X and Y should not be the same as the learning sample,
%   which is the sample that was used to fit the tree T.
%
%   COST = TEST(T,'crossvalidate',X,Y) uses 10-fold cross-validation to
%   compute the cost vector.  X and Y should be the learning sample, which
%   is the sample that was  used to fit the tree T.  The function
%   partitions the sample into 10 subsamples, chosen randomly but with
%   roughly equal size.  For classification trees the subsamples also have
%   roughly the same class proportions.  For each subsample, TEST fits
%   a tree to the remaining data and uses it to predict the subsample.  It
%   pools the information from all subsamples to compute the cost for the
%   whole sample.




if nargin<2, error('stats:treetest:TooFewInputs','Not enough arguments.'); end

if ~ischar(TorCorR) || ~(TorCorR(1)=='t' || TorCorR(1)=='c' || TorCorR(1)=='r')
   error('stats:treetest:InvalidOption',...
         'Second argument must be ''test'', ''crossvalidate'', or ''resubstitution''.');
end
if TorCorR(1)=='t' && nargin<4
   error('stats:treetest:TooFewInputs',...
         'Not enough arguments.  Need X and Y for the test sample.');
elseif TorCorR(1)=='c' && nargin<4
   error('stats:treetest:TooFewInputs',...
         'Not enough arguments.  Need X and Y from the learning sample.');
elseif TorCorR(1)=='r' && nargin>2
    % No param name/value args allowed.  Okay to give X,Y if we recognize
    % X as not a param name
    if nargin>4 || ischar(X)
         error('stats:treetest:TooManyInputs',...
           'Parameter name/value pairs not allowed with the ''resubstitution'' method.');
    end
end
if TorCorR(1)~='r'
   if ~ischar(Y) && numel(Y)~=length(Y)
      error('stats:treetest:BadData','Y must be a vector.');
   else
      if iscell(Y) || isnumeric(Y) || islogical(Y)
         n = length(Y);
      else
         n = size(Y,1);
      end
      if size(X,1)~=n
         error('stats:treetest:InputSizeMismatch',...
               'There must be one Y value for each row of X.');
      end
   end
end

okargs =   {'nsamples' 'treesize'};
defaults = {10         'se'};
[eid,emsg ncv treesize] = dfswitchyard('statgetargs',okargs,defaults,varargin{:});
if ~isempty(emsg)
    error(sprintf('stats:treetest:%s',eid),emsg);
end

if ~isnumeric(ncv) || numel(ncv)~=1 || ncv<2 || ncv~=round(ncv)
   error('stats:treetest:BadNSamples',...
         'Value of ''nsamples'' argument must be an integer 2 or larger.');
end
if ~ischar(treesize) || ~(treesize(1)=='s' || treesize(1)=='m')
   error('stats:treetest:BadTreeSize',...
         'Value of ''treesize'' argument must be ''se'' or ''min''.');
end

% Get complexity parameters for all pruned subtrees
if isempty(Tree.alpha)
   Tree = treeprune(Tree);
end

% Remove missing values
if nargin>=4
   t = any(isnan(X),2);
   
      Yold = Y;
      Y = classname2id(Y,Tree.classname);
      if any(Y==0)
         bad = find(Y==0);
         bad = Yold(bad(1));
         if isnumeric(bad)
            bad = num2str(bad);
         elseif iscell(bad)
            bad = bad{1};
         end
         error('stats:treetest:BadYValue',...
            'At least one Y value (''%s'') is incompatible with the tree.',... 
                       bad);
      end

                   
   t = t | isnan(Y);
   if any(t)
      X(t,:) = [];
      Y(t,:) = [];
   end
end

% Do proper type of testing (error estimation)
switch(TorCorR(1))
 case 't', [cost,secost] = testtree(Tree,X,Y);
 case 'c', [cost,secost] = cvtree(Tree,X,Y,ncv);
 case 'r', [cost,secost] = resubinfo(Tree); treesize = 'm';
end

cost = cost(:);
secost = secost(:);
if nargout>=3
   ntnodes = Tree.ntermnodes(:);
end
if nargout>=4
   bestlevel = selecttree(cost,secost,treesize(1)) - 1;
end

% ---------------------------------------------------------
function [resuberr,seresub] = resubinfo(Tree)
%RESUBINFO Compute error rates for tree using resubstitution error.

% Get complexity parameters for all pruned subtrees
nsub = 1+max(Tree.prunelist);

% Get error rate for each subtree in this sequence
resuberr = zeros(nsub,1);
for j=1:nsub;
   Tj = treeprune(Tree,'level',j-1);
   leaves = Tj.node(Tj.var==0);
   resuberr(j) = sum(Tj.risk(leaves));
end
seresub = zeros(size(resuberr));

% ---------------------------------------------------------------
function [testerr,seerr] = testtree(Tree,X,id)
%TESTTREE Compute error rates for tree using test sample.
%   The id variable is the class id for classification, or the y variable
%   for regression.

% Get pruning sequence and compute fitted values for the whole sequence
nsub = 1 + max(Tree.prunelist);
yfit = treeval(Tree,X,(0:nsub-1));



   nclasses = Tree.nclasses;
   cost = Tree.cost;
   prior = Tree.prior(:);
   if isempty(prior)
      prior = Tree.classcount(1,:)' / Tree.nodesize(1);
   end
   Njtest = histc(id,1:nclasses);
   adjprior = (prior ./ max(eps,Njtest))';


% Compute error statistics

   testerr = zeros(nsub,1);
   seerr = zeros(nsub,1);
   for k=nsub:-1:1;
      % M(i,j) counts class i items classified as class j
      M = accumarray([id,yfit(:,k)], 1, [nclasses,nclasses]);
   
      % Compute loss for this classification
      loss = sum(cost .* M, 2);
      losssq = sum(cost.^2 .* M, 2);
      s2 = losssq  - loss.^2 ./ max(1,Njtest);
   
      testerr(k) = adjprior * loss;
      seerr(k) = sqrt(adjprior.^2 * s2);
   end


% ---------------------------------------------------------------
function [cverr,secverr] = cvtree(Tree,X,id,ncv)
%CVTREE Compute error rates for tree using cross-validation.
%   [CVERR,SECVERR] = CVTREE(TREE,X,ID,NCV)

% Get geometric means of the alpha boundary points
alpha = Tree.alpha;
avgalpha = [sqrt(alpha(1:end-1) .* alpha(2:end)); Inf];

% Loop over cross-validation samples
N = size(X,1);
ntrees = length(avgalpha);
cverr = zeros(ntrees,1);
secverr = zeros(ntrees,1);
cvid = 1 + mod((1:N),ncv);


   % Use a random permutation with fixed category proportions
   idrand = id + rand(size(id));
   [stdid,idx] = sort(idrand);
   cvid = cvid(idx);
   nlevels = size(Tree.cost,1); 
   id = nominal(id,[],1:nlevels);
   coststruct.cost = Tree.cost;
   coststruct.group = 1:nlevels;
   args = {'prior',Tree.prior, 'cost',coststruct, 'prune','on'};


% Get predicted values using cross-validation samples
cvclass = zeros(N,ntrees);
for j=1:ncv
   % Use jth group as a test, train on the others
   
   testrows = find(cvid == j);
   trainrows = find(cvid~=j);
  
   % Get a sequence of pruned trees for the training set
   Tj = multiclasstree(X(trainrows,:),id(trainrows),args{:},'categorical',Tree.catcols,'categorical_group',Tree.catgrp);
   
   % Get classifications based on each subtree that we require
   treesneeded = findsubtree(Tj,avgalpha);
   if isa(Tree,'struct')
    Tj = multiclasstree(Tj);
   end
   [abc abc2]=eval(Tj,X(testrows,:),treesneeded-1);
   abc=str2num(char(abc));
  
   for xyz=1:size(abc2,2)
       
       cvclass(testrows,xyz)=abc(length(testrows)*(xyz-1)+1:length(testrows)*xyz)';
   end
   %cvclass(testrows,:) = double(eval(Tj,X(testrows,:),treesneeded-1));
end

% Compute output statistics based on those predictions

   Nj = Tree.classcount(1,:)';
   prior = Tree.prior;
   if isempty(prior)
      prior = Nj' / N;
   end
   adjprior = (prior ./ max(1,Nj'));
   nclasses = length(prior);
   cost = Tree.cost;
   id = double(id);
   for k=1:ntrees
      M = accumarray([id,cvclass(:,k)], 1, [nclasses,nclasses]);
      loss = sum(cost .* M, 2);
      losssq = sum(cost.^2 .* M, 2);
      s2 = losssq - loss.^2 ./ max(1,Nj);
      cverr(k) = adjprior * loss;
      secverr(k) = sqrt(adjprior.^2 * s2);
   end


% ----------------------------
function k = findsubtree(Tree,alpha0)
%FINDSUBTREE Find subtree corresponding to specified complexity parameters.

adjfactor = 1 + 100*eps;
alpha = Tree.alpha;
k = zeros(size(alpha0));
for j=1:length(alpha0);
   k(j) = sum(alpha <= alpha0(j)*adjfactor);
end

% -----------------------------
function bestj = selecttree(allalpha,sealpha,treesize)
%SELECTTREE Select the best tree from error rates using some criterion.

% Find the smallest tree that gives roughly the minimum error
[minerr,minloc] = min(allalpha);
if isequal(treesize(1),'m')
   cutoff = minerr * (1 + 100*eps);
else
   cutoff = minerr + sealpha(minloc);
end
j = find(allalpha <= cutoff);
bestj = j(end);

% -----------------------------
function idvec = classname2id(idnames,cnames)
%CLASSNAME2ID Create vector of numeric indices from class name array.

idvec = zeros(length(idnames),1);
if isa(idnames,'categorical')
    idvec = double(idnames);
    idnames = char(idnames);
else
    if isnumeric(idnames) || islogical(idnames)
        idnames = cellstr(strjust(num2str(idnames),'left'));
    end
    for j=1:length(cnames)
        idvec(strcmp(cnames(j),idnames)) = j;
    end
end
t = find(idvec==0);
if ~isempty(t)
   txt = idnames(t,:);
   if ischar(txt)
      txt = cellstr(txt);
   end
   idvec(t(cellfun('isempty',txt))) = NaN;
end

