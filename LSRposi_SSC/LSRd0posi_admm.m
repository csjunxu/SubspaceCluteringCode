function C2 = LSRd0posi_admm(Y,affine,Par)

if (nargin < 3)
    % default regularizarion parameters
    Par.lambda = 1;
    Par.rho = 0.1;
    Par.maxIter = 200;
end
if (nargin < 2)
    % default subspaces are linear
    affine = false;
end

% default coefficient and linear system error threshold to stop ADMM
thr = 2*10^-4;

[L, N] = size (Y);


% setting penalty parameters for the ADMM
% mu1 = alpha1 * 1/computeLambda_mat(Y);
% mu2 = alpha2 * 1;

% lambda = computeLambda_mat(Y);

if (~affine)
    %% initialization
    Inv = (Y'*Y+Par.rho/2*eye(N))\eye(N);
    InvW = (2/Par.rho*eye(N)-(2/Par.rho)^2*Y'/(2/Par.rho*(Y*Y')+eye(L))*Y);
    C1 = zeros(N,N);
    Delta = zeros(N,N);
    err1 = 10*thr; err2 = 10*thr;
    i = 1;
    %% ADMM iterations
    while ( err1(i) > thr && i < Par.maxIter )
        %% update A the coefficient matrix
        if N < L
            A = Inv * (Y' * Y + Par.rho/2 * C1 + 0.5 * Delta);
        else
            A =  InvW * (Y' * Y + Par.rho/2 * C1 + 0.5 * Delta);
        end
        A = A - diag(diag(A));
        %% update C the data term matrix
        Q = (Par.rho*A - Delta)/(2*Par.lambda+Par.rho);
        C2  = solver_BCLS_closedForm(Q);
        C2 = C2 - diag(diag(C2));
        
        %% update Deltas the lagrange multiplier matrix
        Delta = Delta + Par.rho * ( C2 - A);
        %% computing errors
        err1(i+1) = errorCoef(C2, A);
        err2(i+1) = errorLinSys(Y, A);
        %%
        C1 = C2;
        i = i + 1;
    end
    fprintf('err1: %2.4f, err2: %2.4f, iter: %3.0f \n',err1(end),err2(end),i);
else
    %% initialization
    Inv = (Y'*Y+Par.rho/2*(eye(N)+ones(N,N)))\eye(N);
    InvW = (2/Par.rho*(eye(N)+ones(N,N))-(2/Par.rho)^2*Y'/(2/Par.rho*(Y*Y')+eye(L))*Y);
    C1 = zeros(N,N);
    Delta = zeros(N,N);
    Delta3 = zeros(1,N);
    err1 = 10*thr; err2 = 10*thr; err3 = 10*thr;
    i = 1;
    %% ADMM iterations
    while ( (err1(i) > thr || err3(i) > thr) && i < Par.maxIter )
        %% update A the coefficient matrix
        if N < L
            A = Inv*(Y'*Y+Par.rho/2*C1+0.5*Delta+Par.rho/2*ones(N,1)*ones(1,N)-ones(N,1)*Delta3);
        else
            A =  InvW*(Y'*Y+Par.rho/2*C1+0.5*Delta+Par.rho/2*ones(N,1)*ones(1,N)-ones(N,1)*Delta3);
        end
        A = A - diag(diag(A));
        %% update C the data term matrix
        Q = (Par.rho*A - Delta)/(2*Par.lambda+Par.rho);
        C2  = solver_BCLS_closedForm(Q);
        C2 = C2 - diag(diag(C2));
        %% updating Lagrange multipliers
        Delta = Delta + Par.rho * (A - C2);
        Delta3 = Delta3 + Par.rho * (ones(1,N)*A - ones(1,N));
        %% computing errors
        err1(i+1) = errorCoef(A,C2);
        err2(i+1) = errorLinSys(Y,A);
        err3(i+1) = errorCoef(ones(1,N)*A,ones(1,N));
        %%
        C1 = C2;
        i = i + 1;
    end
    fprintf('err1: %2.4f, err2: %2.4f, err3: %2.4f, iter: %3.0f \n',err1(end),err2(end),err3(end),i);
end