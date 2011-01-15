function projtest
    
% The goal of these experiments is to see how efficient our algorithm is,
% in particular how the number of iterations required scales with the 
% dimension of the problem.

% Enforce non-negativity?
nn = 1;
    
% These are the various dimensionalities to test    
dims = [2 3 5 10 50 100 500 1000 3000 5000 10000];
    
% These are the various desired sparseness levels
ds = [0.1 0.3 0.5 0.7 0.9];

% These are the various initial sparseness levels
is = [0.1 0.3 0.5 0.7 0.9];

% How many tests for each case?
ntests = 50;

% Initialize matrix to hold results
iters = zeros(length(ds),length(is),length(dims),ntests);

% Go through all different desired sparsenesses
for dsiter = 1:length(ds),
    
    % Go through all different initial sparsenesses
    for isiter = 1:length(is),

	desiredsparseness = ds(dsiter);
	initialsparseness = is(isiter);
	fprintf('[%.1f/%.1f]: ',desiredsparseness,initialsparseness);
	
	% Go through all the different dimensionalities
	for dimiter = 1:length(dims),
    
	    N = dims(dimiter);
	    fprintf('(%d)',N);
    
	    % Take several test cases
	    for testcase = 1:ntests,

		% Take a random vector and project it onto desired sparseness
		x = randn(N,1); x = x/norm(x);
		k1 = sqrt(N)-(sqrt(N)-1)*desiredsparseness; 
		[x,usediters] = projfunc( x, k1, 1, nn );
		
		% Take another random vector and project to initial sparseness
		s = randn(N,1); s = s/norm(s);
		k1 = sqrt(N)-(sqrt(N)-1)*initialsparseness; 
		[s,usediters] = projfunc( s, k1, 1, nn );
		
		% Project s to achieve desired sparseness, save 'usediters'
		k1 = sqrt(N)-(sqrt(N)-1)*desiredsparseness; 		
		[v,usediters] = projfunc( s, k1, 1, nn );
		iters(dsiter,isiter,dimiter,testcase) = usediters;
		
		% If v does not satisfy constraints, then something is wrong!!
		if (abs(sum(abs(v))-k1)>1e-8) | (abs(sum(v.^2)-1)>1e-8),
		    error('L1 or L2 constraint not satisfied!!!!');
		end
		if nn, 
		    if min(v)<0,
			error('Positivity constraint not satisfied!!!');
		    end
		end
		
		% Make sure that s is closer to v than to x!!
		if norm(x-s)<(norm(v-s)-1e-10),
		    error('Not closest point!!!!! Fatal error!!!');
		end
		
	    end
	end
	
	fprintf('\n');
    end
end

% Show average number of iterations as a function of desired sparseness
% and initial sparseness
meaniters = mean(iters,4);
meaniters = mean(meaniters,3);
meaniters

% Note: along the diagonal we really should have zeros, since if the
% initial and the desired sparsenesses match then there is no need for
% even a single iteration! Instead, we have numbers slightly larger than
% one, probably because of roundoff errors.

% Clearly, the 'worst case' is when the desired sparseness is high (0.9)
% and the initial sparseness is very low (0.1). In the paper we plot the
% average number of iterations required for this worst case scenario:
worstcase = reshape(iters(5,1,:,:),[length(dims) ntests]);
meanworstcase = mean(worstcase,2);
maxworstcase = max(worstcase,[],2);
minworstcase = min(worstcase,[],2);
figure(1); 
semilogx(dims,meanworstcase,'k-',dims,maxworstcase,'k:', dims, ...
	 minworstcase,'k:');


% This calculates the maximum overall all trials and conditions
allmaxiters = max(iters(:))

