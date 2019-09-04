function [x, k] = cgp(x0, A, C, b, mit, stol, bbA, bbC)
% Synopsis:
% x0: initial point
% A: Matrix A of the system Ax=b
% C: Preconditioning Matrix can be left or right
% mit: Maximum number of iterations
% stol: residue norm tolerance
% bbA: Black Box that computes the matrix-vector product for A * u
% bbC: Black Box that computes:
%      for left-side preconditioner : ha = C \ ra
%      for right-side preconditioner: ha = C * ra
% x: Estimated solution point
% k: Number of iterations done
%
% Example:
% tic;[x, t] = cgp(x0, S, speye(1), b, 3000, 10^-8, @(Z, o) Z*o, @(Z, o) o);toc
% Elapsed time is 0.550190 seconds.
%
% Reference:
%  Métodos iterativos tipo Krylov para sistema lineales
%  B. Molina y M. Raydan - {{ISBN|908-261-078-X}}
if nargin < 8, error('Not enough input arguments. Try help.'); end;
if isempty(A), error('Input matrix A must not be empty.'); end;
if isempty(C), error('Input preconditioner matrix C must not be empty.'); end;
x = x0;
ha = 0;
hp = 0;
hpp = 0;
ra = 0;
rp = 0;
rpp = 0;
u = 0;
k = 0;

ra = b - bbA(A, x0); % <--- ra = b - A * x0;
while norm(ra, inf) > stol
    ha = bbC(C, ra); % <--- ha = C \ ra;
    k = k + 1;
    if (k == mit), warning('GCP:MAXIT', 'mit reached, no conversion.'); return; end;
    hpp = hp;
    rpp = rp;
    hp = ha;
    rp = ra;
    t = rp' * hp;
    if k == 1
        u = hp;
    else
        u = hp + (t / (rpp' * hpp)) * u;
    end;
    Au = bbA(A, u); % <--- Au = A * u;
    a = t / (u' * Au);
    x = x + a * u;
    ra = rp - a * Au;
end;