function Ds = dictlearn_mb(varargin)
% dictlearn_mb    Learn a dictionary by a minibatch RLS-DLA variant. 
%                 i.e. make a model of the input data X which are available
% through a specified m-file or given as a NxL matrix. 
% The resulting dictionary may be stored in a file, see 'dictfile' option.
% 
%   Use of function:
%   ----------------
%   Ds = dictlearn_mb(varargin)
%     Ds is the dictionary as a struct, Ds.D is the actual dictionary
%     Arguments may be given as pairs: 'option',value, ...
%     (or options may be given in one (or more) struct or cell-array)
% 
%   Options for training data, i.e. how to make the training vectors: 
%   --------------------------------------------------------
%   X        input data, NxL matrix
%            or
%   Xmfile   m-file which return data vectors, i.e. X=feval(Xmfile,Xopt);
%            Options i Xmfile should be processed like in this file
%   Xopt     struct with options for Xmfile, ex: Xopt=struct('L',2000);
%            'L' should be number of training vectors to return
%
%   Options for sparse approximation: 
%   --------------------------------------------------------
%   samet    the sparse approximation method to use. It can be the method
%       argument in sparseapprox.m. Default 'javaormp' 
%       also convex approximation methods: cafixed, convexsets, greedy
%   saopt    the additional options to use in sparseapprox.m as a struct. 
%       Default is: struct('tnz',3, 'verbose',0); 
%   K   the number of atoms in the dictionary, default 100
%
%   Options for minibatch variant of RLS-DLA: 
%   --------------------------------------------------------
%   mb, minibatch   an array with number of batches and number of vectors in each batch.
%       Default: repmat([1000,1],5,1).*[2,10; 2,40; 3,100; 3,200; 2,500];
%       Total number of training vectors is: sum( minibatch(:,1).*minibatch(:,2) )
%   lam0, lambda0   The initial value of the forgetting factor lambda, default 0.996.
%   lam1, lambda1   when lambda should be increased to 1, given as a number relative
%       to the scheduled number of vectors to process, default 0.9.
%   outlim   limit for data outlayers, limit for approximation error relative to
%       mean of the last errors. Default: inf (which is no limit)
%   D    initial dictionary, default: random training vectors
%   A    initial A matrix, default: identity matrix
%
%   Options for mat-files to store dictionary in (or retrieve from): 
%   --------------------------------------------------------
%   dictin, optionfile   name of mat-file where options are stored.
%       This can be a dictionary file.
%   dictfile  file to store Ds in. If not given Ds will not be store
%
%   Options for logging and displaying results: 
%   --------------------------------------------------------
%   ptc, PropertiesToCheck   a cell array with names of dictionary properties to chech
%       at regular intervals during learning given by checkrate. None: {{}}.
%       Default: {{'mtvp','r2avg','traceA','fbA','fbB','betamse','d2min'}}
%       May also use: 'nofr', 'rcondA', 'snr'
%       Note {{ and }} used to define a cell array as a field in a struct!
%       It is good to have 'mtvp' (or 'tvp') as first property (used as x-axis in plots)
%   checkrate   how often properties are reported and dictionary checked, default 1000
%   v, verbose  verbose level given as 0 (only warnings) 1 (normal) 2 (much)
%
%   Examples:
%   ---------
%   X = load('dataXforAR1.mat');
%   opt = struct('K',32, 'samet','mexomp','saopt',struct('tnz',4));
%   Ds = dictlearn_mb('X',X, opt);    
%   figure(1);clf;plot(Ds.ptab(:,1),Ds.ptab(:,2),'b-');xlabel(Ds.ptc{1});ylabel(Ds.ptc{2}); 

%----------------------------------------------------------------------
% Copyright (c) 2013.  Karl Skretting.  All rights reserved.
% University of Stavanger, Signal Processing Group
% Mail:  karl.skretting@uis.no   Homepage:  http://www.ux.uis.no/~karlsk/
% 
% HISTORY:  dd.mm.yyyy
% Ver. 1.0  08.04.2013  Made function (based on texdictlearn.m)
%----------------------------------------------------------------------

mfile = 'DictLearn_MB';

%% defaults, initial values
tstart = tic;
text = [mfile,' started ',datestr(now)];  

Xmfile = '';
Xopt = struct('L',2000);
samet = 'javaORMP';
K = 32;
fac2K = 250;   % only used to set minibatch below
minibatch = repmat([fac2K,1],5,1).*[2,10; 2,40; 3,100; 3,200; 2,500];  
lam0 = 0.99;
lam1 = 0.9;
tvp = 0;    % number of training vectors processed
outnum = 0; 
outlim = inf;
D = [];
A = [];
ptc = {'mtvp','r2avg','traceA','fbA','fbB','betamse','d2min'};
checkrate = 1000;
verbose = 1;
snr = 0;   % used only when X is given

%% get the options
nofOptions = nargin;
optionNumber = 1;
fieldNumber = 1;
while (optionNumber <= nofOptions)
    if isstruct(varargin{optionNumber})
        sOptions = varargin{optionNumber}; 
        sNames = fieldnames(sOptions);
        opName = sNames{fieldNumber};
        opVal = sOptions.(opName);
        % next option is next field or next (pair of) arguments
        fieldNumber = fieldNumber + 1;  % next field
        if (fieldNumber > numel(sNames)) 
            fieldNumber = 1;
            optionNumber = optionNumber + 1;  % next pair of options
        end
    elseif iscell(varargin{optionNumber})
        sOptions = varargin{optionNumber}; 
        opName = sOptions{fieldNumber};
        opVal = sOptions{fieldNumber+1};
        % next option is next pair in cell or next (pair of) arguments
        fieldNumber = fieldNumber + 2;  % next pair in cell
        if (fieldNumber > numel(sOptions)) 
            fieldNumber = 1;
            optionNumber = optionNumber + 1;  % next pair of options
        end
    else
        opName = varargin{optionNumber};
        opVal = varargin{optionNumber+1};
        optionNumber = optionNumber + 2;  % next pair of options
    end
    % interpret opName and opVal
    if  strcmpi(opName,'X') 
        if isnumeric(opVal)
            Xin = opVal;
        else
            error([mfile,': illegal type (class) of value for option ',opName]);
        end
    end
    if  strcmpi(opName,'Xmfile') 
        if ischar(opVal)
            if exist(opVal, 'file')
                Xmfile = opVal;
            else
                error([mfile,': can not find ',opVal]);
            end
        else
            error([mfile,': illegal type (class) of value for option ',opName]);
        end
    end
    if  strcmpi(opName,'Xopt') 
        if isstruct(opVal)
            Xopt = opVal;
        else
            error([mfile,': illegal type (class) of value for option ',opName]);
        end
    end
    %
    if ( strcmpi(opName,'samet') )
        if ischar(opVal)
            samet = opVal;
        else
            error([mfile,': not character value (string) for option ',opName]);
        end
    end 
    if ( strcmpi(opName,'saopt') )
        if isstruct(opVal)
            saopt = opVal;
        else
            error([mfile,': not struct for option ',opName]);
        end
    end 
    if strcmpi(opName,'K')    
        if isnumeric(opVal)
            K = opVal;
        else
            error([mfile,': not numeric for option ',opName]);
        end
    end 
    %
    if ( strcmpi(opName,'minibatch') || strcmpi(opName,'mb') )
        if isnumeric(opVal)
            minibatch = opVal;
        else
            error([mfile,': not numeric for option ',opName]);
        end
    end 
    if ( strcmpi(opName,'lam0') || strcmpi(opName,'lambda0') )
        if isnumeric(opVal)
            lam0 = opVal(1);
        else
            error([mfile,': not numeric for option ',opName]);
        end
    end 
    if ( strcmpi(opName,'lam1') || strcmpi(opName,'lambda1') )
        if isnumeric(opVal)
            lam1 = opVal(1);
        else
            error([mfile,': not numeric for option ',opName]);
        end
    end 
    if ( strcmpi(opName,'r2avg') || strcmpi(opName,'r2average') )
        if isnumeric(opVal)
            r2avg = opVal(1);
        else
            error([mfile,': not numeric for option ',opName]);
        end
    end 
    if ( strcmpi(opName,'outlim') || strcmpi(opName,'outlimit') )
        if isnumeric(opVal)
            outlim = opVal(1);
        else
            error([mfile,': not numeric for option ',opName]);
        end
    end 
    if ( strcmpi(opName,'outnum') || strcmpi(opName,'outnumber') )
        if isnumeric(opVal)
            outnum = opVal(1);
        else
            error([mfile,': not numeric for option ',opName]);
        end
    end 
    if strcmpi(opName,'D') 
        if isnumeric(opVal)
            D = opVal;
        else
            error([mfile,': not numeric for option ',opName]);
        end
    end 
    if strcmpi(opName,'A') 
        if isnumeric(opVal)
            A = opVal;
        else
            error([mfile,': not numeric for option ',opName]);
        end
    end 
    %
    if strcmpi(opName,'dictin') || strcmpi(opName,'optionfile')
        if ischar(opVal)
            dictin = opVal;
            if ((numel(dictin) < 5) || ~strcmpi(dictin((end-3):end),'.mat'))
                dictin = [dictin,'.mat']; %#ok<AGROW>
            end
            if (~exist(dictin, 'file') &&  (numel(strfind(dictin,cat_sep)) == 0))
                dictin = [cat_D,cat_sep,dictin];   %#ok<AGROW>
            end
            if exist(dictin, 'file')
                disp([mfile,': use all fields in ',dictin,' as options.']);
                disp('  !!!  NOTE this is NOT tested yet, it may work though. ');
                % but some few should be ignored
                mfileHERE = mfile;
                tstartHERE = tstart;
                textHERE = text;
                clear text
                load(dictin);       % overwrite variables in the m-file!
                mfile = mfileHERE;  
                tstart = tstartHERE;
                if exist('text','var')
                    text = char(text, textHERE);
                else
                    text = textHERE;
                end
                dictfile = dictin;
                % the more careful code would be like
                %  Ds = load(dictin);  
                %  if isfield(Ds,'imno'); imno = Ds.imno; end;
                %  if isfield(Ds,'imfile'); imfile = Ds.imfile; end;
                %  ....
            else
                error([mfile,': can not find ',dictin]);
            end
        else
            error([mfile,': not char for option ',opName]);
        end
    end
    if strcmpi(opName,'dictfile') || strcmpi(opName,'dictout')
        if ischar(opVal)
            dictfile = opVal;
        else
            error([mfile,': not char for option ',opName]);
        end
    end
    %
    if strcmpi(opName,'tvp') 
        if isnumeric(opVal)
            tvp = opVal(1);
        else
            error([mfile,': not numeric for option ',opName]);
        end
    end 
    if strcmpi(opName,'text') 
        if ischar(opVal)
            text = char(opVal, text);
        else
            error([mfile,': not char for option ',opName]);
        end
    end 
    if ( strcmpi(opName,'PropertiesToCheck') || strcmpi(opName,'ptc') )
        if iscell(opVal)
            ptc = opVal;
        else
            error([mfile,': not cell for option ',opName]);
        end
    end 
    if ( strcmpi(opName,'PropertiesTable') || strcmpi(opName,'ptab') )
        if isnumeric(opVal)
            ptab = opVal;
            pti = size(ptab,1);
        else
            error([mfile,': not numeric for option ',opName]);
        end
    end 
    if ( strcmpi(opName,'checkrate') || strcmpi(opName,'cr') )
        if isnumeric(opVal)
            checkrate = opVal(1);
        else
            error([mfile,': not numeric for option ',opName]);
        end
    end 
    if strcmpi(opName,'verbose') || strcmpi(opName,'v')
        if (islogical(opVal) && opVal); verbose = 1; end;
        if isnumeric(opVal); verbose = opVal(1); end;
    end
end

%% some important variables are checked
K = K(1);
if (~exist('saopt','var') || ~isstruct(saopt) || ~isfield(saopt,'tnz'))
    saopt = struct('tnz',3, 'verbose',0);
end

%% check filename to store Ds into
if exist('dictfile','var') && numel(dictfile)
    if ((numel(dictfile) < 5) || ~strcmpi(dictfile((end-3):end),'.mat'))
        dictfile = [dictfile,'.mat'];
    end
    if exist(dictfile,'file')
        disp([mfile,': dictionary ',dictfile,' exist. It will be overwritten.']);
    end
    if verbose >= 1
        disp([mfile,': samet = ',samet,', dictfile = ''',dictfile,'''.']);
    end
else
    dictfile = '';
end
if exist('Xin','var')   % 
    [N,L] = size(Xin);
else
    if ~(exist(Xmfile,'file') == 2)
        error([mfile,': training data not available.']);
    end
    if ~isfield(Xopt,'L')
        Xopt.L = checkrate;
    end    
    L = Xopt.L;
end

%% initialize variables, including dictionary
StartOfLineDisplayed = false;    % used when (verbose == 1)
if (size(D,1) == N) && (size(D,2) == K) 
    % D is probably ok
else
    if exist('Xin','var')   % 
        D = Xin(:,(L-K+1):L);   % the last ones
    else
        D = feval(Xmfile,Xopt,'L',K);
        if (size(D,2) > K)     % make sure D does not have to many atoms
            D = D(:,1:K);
        end
    end
end
% normalization of dictionary
g = 1./sqrt(sum( D.*D ));
D = D .* repmat(g, N, 1); % normalize dictionary

if (numel(A) == 1)
    B = D/A;   % note A is scalar here
    A = eye(K)*A;
end
if ~((size(A,1) == K) && (size(A,2) == K)) 
    A = eye(K);
    B = D;
end
tvtot = tvp + sum( minibatch(:,1).*minibatch(:,2) );
if ~exist('r2avg','var')
    if exist('Xin','var')   % 
        X = Xin(:,1:(L-K));        % not the last ones (in D) here
    else
        X = feval(Xmfile,Xopt);
    end
    W = sparseapprox(X, D, samet, saopt);
    R = X - D*W;
    r2avg = mean( sqrt(sum(R.*R)) );
end
if ~exist('pti','var');
    pti = 0;      % last written line in ptab
end
if exist('ptab','var') && (size(ptab,2) == numel(ptc)) ;
    ptab = [ptab; zeros(ceil(tvtot/checkrate)+1, numel(ptc))];
elseif (numel(ptc) > 0)
    pti = 0;      % last written line in ptab
    ptab = zeros(ceil(tvtot/checkrate)+1, numel(ptc));
end
ptno = tvp;   % tvp when last line in ptab was written
if (numel(ptc) > 0)   % report some results as learning goes by
    pti = pti + 1;
    % use a function at the end of this m-file, 
    ptab(pti,:) = dictionaryProperties(D, ptc, tvp, r2avg, trace(A), 0, rcond(A), snr);
end

%% Display info
if (verbose > 0)  % verbose
    disp(' ');
    disp([mfile,' started ',datestr(now)]);
end
if (verbose >= 1)  % verbose
    if exist('Xin','var')   % 
        disp([mfile,': training data given as ',int2str(N),'x',int2str(L),' matrix.']);
    else
        disp([mfile,': training data should be made by ',Xmfile,', Xopt.L = ',int2str(Xopt.L)]);
    end
end
if (verbose > 1)  % very verbose
    disp(['  number of dictionary atoms K = ',int2str(K)]);
    disp(['  thus size of dictionary D is ',int2str(size(D,1)),'x',int2str(size(D,2)),'.']);
    disp(['  method for sparse approximation samet = ',samet]);
    disp(['  number of atoms in each sparse representation = ',int2str(saopt.tnz)]);
    disp(['  number of mini-batches to do = ',int2str(sum( minibatch(:,1)))]);
    disp(['  millions training vectors (mtv) to process = ',sprintf('%6.2f',tvtot/1000000)]);
    disp(['  initially lambda = ',num2str(lam0),', it increases to 1 after ',...
           sprintf('%7.3f',lam1*tvtot/1000000),' mtv are processed.']);
    disp(['  outlayer limit, outlim = ',num2str(outlim)]);
    if (numel(ptc) > 0)
        disp(['  for each ',int2str(checkrate),' tv the following dictionary properties are reported:']);
        disp( ['  ',propertyhead(ptc, 10)] );
        disp( ['  ',propertyline(ptab(pti,:), ptc, 10)] );
    end
end

%%  minibatch variant of RLS-DLA
if exist('Xin','var')   % 
    X = Xin;
    %    R = X - D*W;   % R = X initially
    RR = sum(X.*X);    
    sumXX = sum(RR);
    snr = 10*log10(sumXX/sum(RR)); 
else
    X = feval(Xmfile,Xopt);
end
[N,L] = size(X);
xno = 0;
for linje = 1:size(minibatch,1)
    for bno = 1:minibatch(linje,1)
        % if isnan(D(1,1))
        %     error(['D(1,1) is NaN, tvp=',int2str(tvp),', linje=',int2str(linje),', bno=',int2str(bno)]);
        % end
        batchsize = minibatch(linje,2);
        if exist('Xin','var')   % 
            if rand(1) < 0.25
                idx = rem(xno+(1:batchsize)-1,L)+1;
                xno = idx(end);
            else
                idx = ceil(rand(1,batchsize)*L);   % just random vectors
            end
            Xbatch = X(:,idx);
        else
            if (xno+batchsize) > size(X,2)        % get more training vectors
                X = feval(Xmfile,Xopt);
                xno = 0;
            end
            Xbatch = X(:,xno+(1:batchsize));
            xno = xno+batchsize;
        end
        %
        W = sparseapprox(Xbatch, D, samet, saopt);
        temp = nnz(isnan(W));
        if temp
            disp(' ');
            disp(['  (479) ',int2str(temp),' elements of W are NaN.  Should stop now.']);
            pause(60)
        end
        temp = nnz(isnan(D));
        if temp
            disp(' ');
            disp(['  (485) ',int2str(temp),' elements of D are NaN.  Should stop now.']);
            pause(60)
        end
        R = Xbatch - D*W;
        if exist('Xin','var')   % 
            RR(idx) = sum(R.*R);
            snr = 10*log10(sumXX/sum(RR)); 
        end
        r2 = sqrt(sum(R.*R));   % 2-norm of approximation error
        if (outlim < inf)
            I = (r2 <= (outlim*r2avg));
            Xbatch = Xbatch(:,I);
            W = W(:,I);
            r2 = r2(I);
            batchsize = numel(r2);
            outnum = outnum + (numel(I)-nnz(I));
        end
        r2avg = 0.8*r2avg + 0.2*mean(r2);
        %
        lam = lambdafun(tvp, 'Q', lam1*tvtot, lam0, 1).^batchsize;
        A = lam*A + full(W*W');
        B = lam*B + full(Xbatch*W');
        try
            D = B/A;    % it may happen that A is ill conditioned
            temp = nnz(isnan(D));
            if temp
                disp(' ');
                disp(['  (499) ',int2str(temp),' elements of D are NaN.  Should stop now.']);
                pause(60)
            end
            % then a warning is displayed and no error thrown
            % the "normal" place to adjust dictionary and A is ca 80 lines below.
        catch ME
            if StartOfLineDisplayed
                disp('  '); 
                StartOfLineDisplayed = false;
            end
            disp(['  509: ',ME.message]);
            if strcmpi(samet((end-1):end),'MP') 
                % here we replace 'unused' atoms by atoms from training set
                d = diag(A);
                counter = 0;
                while (min(d) < 0.1*mean(d))  && (counter < min(20, size(Xbatch,2)))
                    counter = counter + 1;
                    [~, atom_no] = min(d);
                    D(:,atom_no) = Xbatch(:,counter); 
                    A(atom_no,:) = 0;
                    A(:,atom_no) = 0;
                    A(atom_no,atom_no) = 0.5*mean(d);
                    d = diag(A);
                end
                B = D*A;
                disp(['  replace ',int2str(counter),' atoms in D by training vectors.']);
            else % start with 'initial' values again
                if exist('Xin','var')   %
                    D = Xin(:,ceil(L*rand(1,K))); % K random numbers
                else
                    D = feval(Xmfile,Xopt,'L',K);
                    if size(D,2) > K
                        D = D(:,1:K);
                    end
                end
                g = 1./sqrt(sum( D.*D ));
                D = D .* repmat(g, N, 1); % normalize dictionary
                D = D + 0.1*randn(N,K);   % add some noise
                A = eye(K);
                B = D;
                disp('  replace the whole dictionary D by training vectors.');
            end
        end
        % normalization of dictionary
        g = 1./sqrt(sum( D.*D ));
        D = D .* repmat(g, N, 1); % normalize dictionary, not A and B though
        tvp = tvp+batchsize;
        %
        if (numel(ptc) > 0) % report results (ptc: properties to check)
            if (tvp >= (ptno + checkrate))
                % 21.feb.2012: moved the report part after the check part
                pti = min(pti+1, size(ptab,1));
                ptno = tvp;   % tvp when last line in ptab was written
                %
                % the matrix A is also checked at checkrate
                d = diag(A);
                md = mean(d);
                if (min(d) < 0);   % should not happen, numerical error??
                    if StartOfLineDisplayed
                        disp('  '); 
                        StartOfLineDisplayed = false;
                    end
                    disp('  negative element in diagonal of A is an error, but will be corrected here.');
                    A = eye(K);   
                    B = D;
                    nof_replaced = K;
                    rcond_A = 1;
                elseif (min(d) < 0.01*md)
                    % small diagonal value mean an atom not much used
                    % adjust this in different ways depending on the used vector seletion method.
                    % display message only if (verbose > 1)
                    if (strcmpi(samet, 'cafixed') || strcmpi(samet, 'convexsets'))
                        nof_replaced = 0;
                        while (min(d) < 0.01*md)  && (nof_replaced < min(K/4, size(Xbatch,2)))
                            nof_replaced = nof_replaced + 1;
                            [~, atom_no] = min(d);
                            neighbor_rows = (sum(sety == atom_no,2) > 0);
                            neighbors = unique(sety( repmat(neighbor_rows,1,size(sety,2)) ));
                            neighbors = setdiff(neighbors, atom_no);
                            if (numel(neighbors) > 5) % use 5 of the neighbors
                                temp = randperm(numel(neighbors));
                                neighbors = neighbors(temp(1:5));
                            end
                            D(:,atom_no) = mean(D(:,neighbors),2);  
                            D(:,atom_no) = (3*D(:,atom_no) +  1*Xbatch(:,nof_replaced))/4;  
                            A(atom_no,:) = 0;
                            A(:,atom_no) = 0;
                            A(atom_no,atom_no) = 0.5*mean(d(neighbors));
                            d = diag(A);
                            % disp([' Replace atom ',int2str(atom_no),' by mean of neighbors']);
                            % disp(int2str(neighbors));
                            % disp([' and 1/4 from training vector ',int2str(nof_replaced),' in Xbatch.']);
                            % temp = nnz(isnan(D));
                            % if temp
                            %     disp(' ');
                            %     disp(['  (740) ',int2str(temp),' elements of D are NaN.  Should stop now.']);
                            %     pause(60)
                            % end
                        end
                        if (verbose > 1) && numel(find(strcmpi(ptc, 'nofr')))
                            fprintf('  %9.3f    replace %i atoms in D by convex combinations of neighbors amd tv.\n',tvp/1e6, nof_replaced);
                        end
                        %
                        temp = nnz(isnan(D));
                        if temp
                            disp(' ');
                            disp(['  (605) ',int2str(temp),' elements of D are NaN.  Should stop now.']);
                            pause(60)
                        end
                        %
                    else
                        I = ((100*d) < md);
                        nof_replaced = nnz(I);
                        if (nof_replaced > 0)
                            if (verbose > 1) && numel(find(strcmpi(ptc, 'nofr')))
                                % fprintf('\n  %9.3f    replace %i atoms with training vectors.\n',tvp/1e6,nof_replaced);
                                fprintf('\n  %9.3f    replace %i atoms with other atoms.\n',tvp/1e6,nof_replaced);
                            end
                            [~,temp] = sort(d,'descend');
                            % D(:,I) = X(:, ceil(rand(1,nof_replaced)*size(X,2))); %#ok<AGROW>
                            D(:,I) = D(:,temp(1:nof_replaced)) + 0.1*randn(N,nof_replaced); 
                            A(:,I) = zeros(K,nof_replaced);
                            A(I,:) = zeros(nof_replaced,K);
                            A(I,I) = (0.5*md)*eye(nof_replaced);
                            if ~(strcmpi(samet, 'greedy'))
                                % here we have:  ~strcmpi(samet, 'cafixed'), and
                                % normalize the dictionary D
                                g = 1./sqrt(sum( D.*D ));
                                D = D .* repmat(g, N, 1); % normalize dictionary, not A and B though
                            end
                        end
                    end
                    %  only when min(d) is small rcond(A) is checked
                    rcond_A = rcond(A);
                    if rcond_A < 1e-8
                        A = A +(0.0001*md)*eye(K);
                        if (verbose > 1) && (numel( find(strcmpi(ptc, 'rcondA')) ) == 0)
                            fprintf('  %9.3f    adjust A to reduce rcond(A) from %g to %g.\n',tvp/1e6, rcond_A, rcond(A));
                        end
                    end
                    %
                    B = D*A;
                    %
                else  % min(d) is large enough 
                    nof_replaced = 0;
                    rcond_A = rcond(A);
                end
                %
                % report to property table (logg),
                ptab(pti,:) = dictionaryProperties(D, ptc, tvp, r2avg, trace(A), ...
                    nof_replaced, rcond_A, snr);
                %
                if (verbose == 1)
                    % if you want to pay extra attention during learning it may
                    % be helful to see the value of rcond(A), this is often
                    % getting very close to zero (or infinity) when learning
                    % goes off in the wrong direction.
%                     if ~StartOfLineDisplayed
%                         fprintf('  Processing up to %7.3f mtv, mtvp = %7.3f (%7.3f)',...
%                             tvtot/1e6, tvp/1e6,log10(rcond_A));
%                         StartOfLineDisplayed = true;
%                     else
%                         fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%7.3f (%7.3f)',...
%                             tvp/1e6,log10(rcond_A));
%                     end
                    % but under normal condition it is just helpful to se number 
                    % of iterations  processed, to assure that learning is
                    % progressing as expected
                    if ~StartOfLineDisplayed
                        fprintf('  Processing up to %7.3f mtv, mtvp = %7.3f ',...
                            tvtot/1e6, tvp/1e6 );
                        StartOfLineDisplayed = true;
                    else
                        fprintf('\b\b\b\b\b\b\b\b%7.3f ', tvp/1e6);
                    end
                elseif (verbose > 1)
                    if (rem(pti,10)==0); disp( ['   ',propertyhead(ptc, 10)] ); end;
                    disp( ['  ',propertyline(ptab(pti,:), ptc, 10)] );
                end
                %
            end
        else  % nothing to report (or check! obs!) 
            % used when fast execution is wanted
        end
        %
    end
end
if exist('ptab','var') && (pti < size(ptab,1)) 
    if (tvp > ptno)  
        pti = pti + 1;
        % ptno = tvp;    does not need this more
        ptab(pti,:) = dictionaryProperties(D, ptc, tvp, r2avg, trace(A), 0, rcond(A), snr);
    end
    ptab = ptab(1:pti,:);
end
if ~exist('ptab','var')
    ptab = [];
end

%% set output argument
tu = toc(tstart);
tuh = floor(tu/3600);
tum = floor(tu/60) - tuh*60;
tus = ceil(tu) - tum*60 - tuh*3600;
t1 = sprintf('%s finished %s, time used is %i:%02i:%02i (hh:mm:ss), %6.3f ms/tvp.', ...
         mfile, datestr(now), tuh, tum, tus, tu*1000/tvp );
text = char(text, t1);
Ds = struct( 'dictfile', dictfile, 'mfile', mfile, ...
             'Xmfile', Xmfile, 'Xopt', Xopt, ...
             'samet', samet, 'saopt', saopt, ...
             'lam0', lam0, 'lam1', lam1, ...
             'N', N, 'K', K, 'L',L, 'D', D, ...
             'r2avg', r2avg, ...
             'ptc', {ptc}, 'ptab', ptab, 'checkrate', checkrate, ...
             'minibatch', minibatch, 'tvp', tvp, ...
             'timeused', tu, ...
             'text', text );
%
if numel(dictfile)
    save(dictfile,'-struct','Ds');
end


%% finish and return

if (verbose > 0) 
    disp(' ');
    disp(t1);
    if (verbose > 1)  % very verbose
        disp(['  number of training vectors processed tvp = ',int2str(tvp)]);
        disp(['  number of training vectors rejected = ',int2str(outnum)]);
    end
    disp(' ');
end

return


function v = dictionaryProperties(D, propnames, tvp, r2avg, traceA, nofr, rcondA, snr)
v = zeros(1,numel(propnames));
m = find(strcmpi(propnames, 'tvp'));
if (numel(m)==1); v(m) = tvp; end;
m = find(strcmpi(propnames, 'mtvp'));
if (numel(m)==1); v(m) = tvp/1e6; end;
m = find(strcmpi(propnames, 'r2avg'));
if (numel(m)==1); v(m) = r2avg; end;
m = find(strcmpi(propnames, 'traceA'));
if (numel(m)==1); v(m) = traceA; end;
m = find(strcmpi(propnames, 'nofr'));
if (numel(m)==1); v(m) = nofr; end;
m = find(strcmpi(propnames, 'rcondA'));
if (numel(m)==1); v(m) = rcondA; end;
m = find(strcmpi(propnames, 'snr'));
if (numel(m)==1); v(m) = snr; end;
if (sum(sum(isnan(D))) > 0)
    disp('  DictLearn_MB.dictionaryProperties (894): Warning D is NaN.');
    return;
end
p = dictprop(dictnormalize(D), 0);
m = find(strcmpi(propnames, 'mu'));
if (numel(m)==1); v(m) = p.mu; end;
m = find(strcmpi(propnames, 'mumin'));
if (numel(m)==1); v(m) = p.mumin; end;
m = find(strcmpi(propnames, 'muavg'));
if (numel(m)==1); v(m) = p.muavg; end;
m = find(strcmpi(propnames, 'mumse'));
if (numel(m)==1); v(m) = p.mumse; end;
m = find(strcmpi(propnames, 'fbA'));
if (numel(m)==1); v(m) = p.A; end;
m = find(strcmpi(propnames, 'fbB'));
if (numel(m)==1); v(m) = p.B; end;
m = find(strcmpi(propnames, 'betamin'));
if (numel(m)==1); v(m) = p.betamin; end;
m = find(strcmpi(propnames, 'betaavg'));
if (numel(m)==1); v(m) = p.betaavg; end;
m = find(strcmpi(propnames, 'betamse'));
if (numel(m)==1); v(m) = p.betamse; end;
%
K = size(D,2);
if sum(strcmpi(propnames, 'd1min') | strcmpi(propnames, 'd1avg')) > 0
    % 1-norm distance between atoms 
    A = zeros(K);
    for k=1:K
        A(k,:) = sum(abs( repmat( D(:,k), 1, K) - D ));
    end
    m = find(strcmpi(propnames, 'd1avg'));
    if (numel(m)==1); v(m) = sum(A(:))/(K*(K-1)); end;
    A = A + max(A(:))*eye(K);
    m = find(strcmpi(propnames, 'd1min'));
    if (numel(m)==1); v(m) = min(A(:)); end;
end
if sum(strcmpi(propnames, 'd2min') | strcmpi(propnames, 'd2avg')) > 0
    % 2-norm distance between atoms 
    A = zeros(K);
    for k=1:K
        A(k,:) = sqrt(sum( ( repmat( D(:,k), 1, K) - D ).^2 ));
    end
    m = find(strcmpi(propnames, 'd2avg'));
    if (numel(m)==1); v(m) = sum(A(:))/(K*(K-1)); end;
    A = A + max(A(:))*eye(K);
    m = find(strcmpi(propnames, 'd2min'));
    if (numel(m)==1); v(m) = min(A(:)); end;
end

return

function s = propertyhead(propnames, len)
s = blanks(len*numel(propnames));
ii = 1;
for i=1:numel(propnames);
    i2 = floor((len-length(propnames{i}))/2);
    i1 = len-length(propnames{i})-i2;
    s(ii:(ii+len-1)) = [blanks(i1),propnames{i},blanks(i2)];
    ii = ii+len;
end
return

function s = propertyline(v, propnames, len)
mudes = '.4f ';
betades = '.2f ';
s = blanks(len*numel(propnames));
ls = int2str(len-1);
ls3 = int2str(len-3);
m = find(strcmpi(propnames, 'tvp'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,'i '],v(m));
end
m = find(strcmpi(propnames, 'snr'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,betades],v(m));
end
m = find(strcmpi(propnames, 'nofr'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls3,'i   '],v(m));
end
m = find(strcmpi(propnames, 'mtvp'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,'.3f '],v(m));
end
m = find(strcmpi(propnames, 'r2avg'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,'.3f '],v(m));
end
m = find(strcmpi(propnames, 'traceA'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,'.2e '],v(m));
end
m = find(strcmpi(propnames, 'rcondA'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,'.2e '],v(m));
end
m = find(strcmpi(propnames, 'mu'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,mudes],v(m));
end
m = find(strcmpi(propnames, 'mumin'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,mudes],v(m));
end
m = find(strcmpi(propnames, 'muavg'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,mudes],v(m));
end
m = find(strcmpi(propnames, 'mumse'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,mudes],v(m));
end
m = find(strcmpi(propnames, 'fbA'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,'.4f '],v(m));
end
m = find(strcmpi(propnames, 'fbB'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,'.1f '],v(m));
end
m = find(strcmpi(propnames, 'betamin'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,betades],v(m));
end
m = find(strcmpi(propnames, 'betaavg'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,betades],v(m));
end
m = find(strcmpi(propnames, 'betamse'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,betades],v(m));
end
m = find(strcmpi(propnames, 'd1min'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,'.3f '],v(m));
end
m = find(strcmpi(propnames, 'd1avg'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,'.3f '],v(m));
end
m = find(strcmpi(propnames, 'd2min'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,'.3f '],v(m));
end
m = find(strcmpi(propnames, 'd2avg'));
if (numel(m)==1)
    ii = 1+(m-1)*len;
    s(ii:(ii+len-1)) = sprintf(['%',ls,'.3f '],v(m));
end
return
