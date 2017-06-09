function C2 = LSRd0posi_admmOutlier(Y,affine,Par)

if (nargin < 3)
    % default regularizarion parameters
    Par.alpha = 20;
    Par.rho = 0.1;
    Par.maxIter = 200;
end

if (nargin < 2)
    % default subspaces are linear
    affine = false; 
end

    % default coefficient error threshold to stop ALM
    % default linear system error threshold to stop ALM
    thr = 2*10^-4; 

if (length(alpha) == 1)
    alpha1 = Par.alpha(1);
    alpha2 = Par.alpha(1);
    alpha3 = Par.alpha(1);
elseif (length(alpha) == 2)
    alpha1 = Par.alpha(1);
    alpha2 = Par.alpha(2);
    alpha3 = Par.alpha(2);
elseif (length(alpha) == 3)
    alpha1 = Par.alpha(1);
    alpha2 = Par.alpha(2);
    alpha3 = Par.alpha(3);
end


[D,N] = size(Y);

gamma = alpha3 / norm(Y,1);
P = [Y eye(D)/gamma];

% setting penalty parameters for the ADMM
mu1 = alpha1 * 1/computeLambda_mat(Y,P);
mu2 = alpha2 * 1;

W = ones(N+D,N);

if (~affine)
    % initialization
%     A = inv(mu1*(P'*P)+mu2*eye(N+D));
    Inv = (mu1*(P'*P)+mu2*eye(N+D))\eye(N+D);
    InvW = (1/mu2*eye(N+D) - (1/mu2)^2*P'/(1/mu2*(P*P')+ eye(D))*P);
    C1 = zeros(N+D,N);
    Lambda1 = zeros(D,N);
    Lambda2 = zeros(N+D,N);
    err1 = 10*thr; err2 = 10*thr;
    i = 1;
    % ADMM iterations
    while ( i < Par.maxIter )%(err1(i) > thr1 || err2(i) > thr2) &&
        % updating Z
        Z = A * (mu1*P'*(Y+Lambda1/mu1)+mu2*(C1-Lambda2/mu2));
        Z(1:N,:) = Z(1:N,:) - diag(diag(Z(1:N,:)));
        % updating C
        C2 = max(0,(abs(Z+Lambda2/mu2) - 1/mu2*W)) .* sign(Z+Lambda2/mu2);
        C2(1:N,:) = C2(1:N,:) - diag(diag(C2(1:N,:)));
        % updating Lagrange multipliers
        Lambda1 = Lambda1 + mu1 * (Y - P * Z);
        Lambda2 = Lambda2 + mu2 * (Z - C2);
        % computing errors
        err1(i+1) = errorCoef(Z,C2);
        err2(i+1) = errorLinSys(P,Z);
        %
        C1 = C2;
        i = i + 1;
    end
    fprintf('err1: %2.4f, err2: %2.4f, iter: %3.0f \n',err1(end),err2(end),i);
else
    % initialization
    delta = [ones(N,1);zeros(D,1)];
    A = inv(mu1*(P'*P)+mu2*eye(N+D)+mu2*(delta*delta'));
    C1 = zeros(N+D,N);
    Lambda1 = zeros(D,N);
    Lambda2 = zeros(N+D,N);
    lambda3 = zeros(1,N);
    err1 = 10*thr1; err2 = 10*thr2; err3 = 10*thr1;
    i = 1;
    % ADMM iterations
    while ( (err1(i) > thr1 || err2(i) > thr2 || err3(i) > thr1) && i < maxIter )
        % updating Z
        Z = A * (mu1*P'*(Y+Lambda1/mu1)+mu2*(C1-Lambda2/mu2)+mu2*delta*(ones(1,N)-lambda3/mu2));
        Z(1:N,:) = Z(1:N,:) - diag(diag(Z(1:N,:)));
        % updating C
        C2 = max(0,(abs(Z+Lambda2/mu2) - 1/mu2*ones(N+D,N))) .* sign(Z+Lambda2/mu2);
        C2(1:N,:) = C2(1:N,:) - diag(diag(C2(1:N,:)));
        % updating Lagrange multipliers
        Lambda1 = Lambda1 + mu1 * (Y - P * Z);
        Lambda2 = Lambda2 + mu2 * (Z - C2);
        lambda3 = lambda3 + mu2 * (delta'*Z - ones(1,N));
        % computing errors
        err1(i+1) = errorCoef(Z,C2);
        err2(i+1) = errorLinSys(P,Z);
        err3(i+1) = errorCoef(delta'*Z,ones(1,N));
        %
        C1 = C2;
        i = i + 1;
    end
    fprintf('err1: %2.4f, err2: %2.4f, err3: %2.4f, iter: %3.0f \n',err1(end),err2(end),err3(end),i);
end