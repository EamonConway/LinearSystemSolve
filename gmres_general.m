function dx = gmres_general(mvp,b,x0,m,psolve,tol)
% mvp - Matrix Vector Product - is a function handle to the implementation of the matrix vector product. 
% b - Right hand side vector.
% x0 - Initial solution guess, 
% m - restart value.
% psolve - PRECONDITIONER SOLVE - is a function handle to P^(-1)b = x
% tol - tolerance of solve.


% Initial implemenation of multigrid. 
% global mg_Jacobian_L
% [~,q] = size(mg_Jacobian_L);
numRow = length(x0);
dx        = x0;
residual = inf;
U = zeros(numRow,m);
H = zeros(m,m);
c = zeros(m,1);
s = c;

%%%% Arnoldi Iterations %%%%

%% Run
    func_evals = 0;
while (residual>tol)
% fprintf('Initialising Krylov subspace \n')
x0     = dx;
r0     = b - mvp(x0);
func_evals = func_evals+1;
beta   = norm(r0);
Q(:,1) = r0/beta;

gammasol = zeros(m+1,1);
gammasol(1)= beta;  

for j  = 1:m
    U(:,j) = psolve(Q(:,j));
    w      = mvp(U(:,j));
    func_evals = func_evals+1;
    
    for i  = 1:j
        
        H(i,j) = w'*Q(:,i);       
        w      = w - H(i,j)*Q(:,i);
      
    end
    
    H(j+1,j)   = norm(w);
    
    if (H(j+1,j)<eps)
      
        break
        
    end
    
Q(:,j+1) = w/H(j+1,j);

  for i = 1:j-1
            H(i:i+1,j)= [c(i), s(i);-s(i), c(i)]*H(i:i+1,j);
  end
        % solve for c and s
        % use for givens
        c(j) = H(j,j)/sqrt(H(j,j)^2 + H(j+1,j)^2);
        s(j) = H(j+1,j)/sqrt(H(j,j)^2 + H(j+1,j)^2);
        %perform givens rotation on H and g
        H(j:j+1,j)= [c(j), s(j);-s(j), c(j)]*H(j:j+1,j);
        gammasol(j:j+1) = [c(j)*gammasol(j); - s(j)*gammasol(j)];

        residual = norm(gammasol(j+1,1));%%uses givens rotations to find residual

if (residual<=tol) 
        break    % Cause Im Happyyyyy, clap along if you feeeeeel
end

end
yk   = H(1:j,1:j)\gammasol(1:j,1);

dxadd = U(:,1:j)*yk;
dx    = x0 + dxadd;

end

fprintf('GMres requires %d function evaluations \n',func_evals)

end        