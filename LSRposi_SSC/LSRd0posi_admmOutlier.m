function C2 = LSRd0posi_admmOutlier(Y,affine,Par)

if (nargin < 3)
    % default regularizarion parameters
    Par.rho = 20;
    Par.maxIter = 200;
    % default coefficient error threshold to stop ALM
    % default linear system error threshold to stop ALM
    Par.thr = 2*10^-4;
end

if (nargin < 2)
    % default subspaces are linear
    affine = false;
end

if (length(rho) == 1)
    rho1 = Par.rho(1);
    rho2 = Par.rho(1);
    rho3 = Par.rho(1);
elseif (length(rho) == 2)
    rho1 = Par.rho(1);
    rho2 = Par.rho(2);
    rho3 = Par.rho(2);
elseif (length(rho) == 3)
    rho1 = Par.rho(1);
    rho2 = Par.rho(2);
    rho3 = Par.rho(3);
end



[D,N] = size(Y);

gamma = rho3 / norm(Y,1);
P = [Y eye(D)/gamma];
[D,L] = size(P);
% setting penalty parameters for the ADMM
mu1 = rho1 * 1/computeLambda_mat(Y,P);
mu2 = rho2 * 1;

if (~affine)
    %% initialization
    Inv = (mu1*(P'*P)+mu2/2*eye(N+D))\eye(N+D);
    InvW = (2/Par.rho*eye(N+D) - (2/Par.rho)^2*P'/(2/Par.rho*(P*P')+ eye(D))*P);
    C1 = zeros(N+D,N);
    Delta1 = zeros(D,N);
    Delta2 = zeros(N+D,N);
    err1 = 10*Par.thr; err2 = 10*Par.thr;
    i = 1;
    %% ADMM iterations
    while ( err1(i) > thr || err2(i) > thr && i < Par.maxIter )
        %% update A the coefficient matrix
        if L < D
            A = Inv*(P'*(Y+Par.rho/2*Delta1)+Par.rho/2*C1+0.5 * Delta2);
        else
            A =  InvW*(P'*(Y+Par.rho/2*Delta1)+Par.rho/2*C1+0.5 * Delta2);
        end
        A(1:N,:) = A(1:N,:) - diag(diag(A(1:N,:)));
        %% update C the data term matrix
        Q = (Par.rho*A - Delta2)/(2*Par.lambda+Par.rho);
        C2  = solver_BCLS_closedForm(Q);
        C2(1:N,:) = C2(1:N,:) - diag(diag(C2(1:N,:)));
        %% updating Lagrange multipliers
        Delta1 = Delta1 + mu1 * (Y - P * A);
        Delta2 = Delta2 + rho2 * (C2 - A);
        %% computing errors
        err1(i+1) = errorCoef(A, C2);
        err2(i+1) = errorLinSys(P, A);
        %%
        C1 = C2;
        i = i + 1;
        fprintf('err1: %2.4f, err2: %2.4f, iter: %3.0f \n',err1(end),err2(end),i);
    end
    fprintf('err1: %2.4f, err2: %2.4f, iter: %3.0f \n',err1(end),err2(end),i);
else
    %% initialization
    delta = [ones(N,1);zeros(D,1)];
    Inv = (P'*P+Par.rho/2*eye(N+D)+Par.rho/2*(delta*delta'))\eye(N+D);
    InvW = (2/Par.rho*(eye(N+D)+(delta*delta'))-(2/Par.rho)^2*P'/(2/Par.rho*(P*P')+ eye(D))*P);
    C1 = zeros(N+D,N);
    Delta1 = zeros(D,N);
    Delta2 = zeros(N+D,N);
    Delta3 = zeros(1,N);
    err1 = 10*Par.thr; err2 = 10*Par.thr; err3 = 10*Par.thr;
    i = 1;
    %% ADMM iterations
    while ( (err1(i) > thr1 || err2(i) > thr2 || err3(i) > thr1) && i < maxIter )
        %% updating Z
        if L < D
            A = Inv*(P'*(Y+Par.rho/2*Delta1)+Par.rho/2*C1+0.5*Delta2+Par.rho/2*delta*ones(1,N)-0.5*delta*Delta3);
        else
            A =  InvW*(P'*(Y+Par.rho/2*Delta1)+Par.rho/2*C1+0.5*Delta2+Par.rho/2*delta*ones(1,N)-0.5*delta*Delta3);
        end
        A(1:N,:) = A(1:N,:) - diag(diag(A(1:N,:)));
        %% update C the data term matrix
        Q = (Par.rho*A - Delta2)/(2*Par.lambda+Par.rho);
        C2  = solver_BCLS_closedForm(Q);
        C2(1:N,:) = C2(1:N,:) - diag(diag(C2(1:N,:)));
        %% updating Lagrange multipliers
        Delta1 = Delta1 + Par.rho * (Y - P * A);
        Delta2 = Delta2 + Par.rho * (C2 - A);
        Delta3 = Delta3 + Par.rho * (delta'*A - ones(1,N));
        %% computing errors
        err1(i+1) = errorCoef(A, C2);
        err2(i+1) = errorLinSys(P, A);
        err3(i+1) = errorCoef(delta'*A, ones(1,N));
        %%
        C1 = C2;
        i = i + 1;
    end
    fprintf('err1: %2.4f, err2: %2.4f, err3: %2.4f, iter: %3.0f \n',err1(end),err2(end),err3(end),i);
end