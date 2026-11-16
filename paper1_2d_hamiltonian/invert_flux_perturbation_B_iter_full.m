function [dg, Gext, p, info] = invert_flux_perturbation_B_iter_full( ...
    obs, g_base, Ni, Ne, lambda, sigma_r, sigma_B, maxit)
% Iterative Gauss-Newton inversion with endpoint |B| constraints and the
% full internal perturbation vector dg_1..dg_Ni (including the dipole).
if nargin < 8 || isempty(maxit), maxit = 20; end
Kint = Ni; Kext = Ne; K = Kint + Kext;

% --- path rows -------------------------------------------------------------
Apath = []; bpath = [];
for j = 1:numel(obs)
    r_obs = obs(j).r_obs; theta_obs = obs(j).theta;
    M = numel(r_obs) - 1;
    r_a = r_obs(1); th_a = theta_obs(1);
    Xi0_a = 0;
    for n = 1:Ni
        [Phi,~,~] = flux_basis(r_a, th_a, n);
        Xi0_a = Xi0_a + g_base(n)*Phi;
    end
    Aj = zeros(M,K); bj = zeros(M,1);
    for i = 2:numel(r_obs)
        r_i = r_obs(i); th_i = theta_obs(i);
        Xi0_i = 0;
        for n = 1:Ni
            [Phi,~,~] = flux_basis(r_i, th_i, n);
            Xi0_i = Xi0_i + g_base(n)*Phi;
        end
        bj(i-1) = -(Xi0_i - Xi0_a);
        for k = 1:Kint
            n = k;
            [Phi_i,~,~] = flux_basis(r_i, th_i, n);
            [Phi_a,~,~] = flux_basis(r_a, th_a, n);
            Aj(i-1,k) = Phi_i - Phi_a;
        end
        for n = 1:Kext
            [Phi_i,~,~] = flux_basis_ext(r_i, th_i, n);
            [Phi_a,~,~] = flux_basis_ext(r_a, th_a, n);
            Aj(i-1,Kint+n) = Phi_i - Phi_a;
        end
    end
    Apath = [Apath; Aj]; %#ok<AGROW>
    bpath = [bpath; bj]; %#ok<AGROW>
end

% --- initial guess: path-only full inversion -------------------------------
[dg0, G0] = invert_flux_perturbation_full(obs, g_base, Ni, Ne, lambda);
p = [dg0; G0];

for it = 1:maxit
    dg = p(1:Kint); Gext = p(Kint+1:end);
    gint = g_base + dg;

    H = []; rhs = [];
    for j = 1:numel(obs)
        pts = [obs(j).r_obs(1), obs(j).theta(1), obs(j).B_start_obs; ...
               obs(j).r_obs(end), obs(j).theta(end), obs(j).B_end_obs];
        for ip = 1:2
            r_p = pts(ip,1); th_p = pts(ip,2); B_obs = pts(ip,3);
            [~, Br0, Bt0, B0] = flux_function_ext(r_p, th_p, gint, Gext);
            bhat = [Br0; Bt0]/B0;
            row = zeros(1,K);
            for k = 1:Kint
                n = k;
                [~, dPr, dPt] = flux_basis(r_p, th_p, n);
                Br = dPt/(r_p^2*sin(th_p));
                Bt = -dPr/(r_p*sin(th_p));
                row(k) = bhat(1)*Br + bhat(2)*Bt;
            end
            for n = 1:Kext
                [~, dPr, dPt] = flux_basis_ext(r_p, th_p, n);
                Br = dPt/(r_p^2*sin(th_p));
                Bt = -dPr/(r_p*sin(th_p));
                row(Kint+n) = bhat(1)*Br + bhat(2)*Bt;
            end
            H = [H; row/sigma_B]; %#ok<AGROW>
            rhs = [rhs; (B_obs-B0)/sigma_B]; %#ok<AGROW>
        end
    end

    A = [Apath; H];
    r = [bpath - Apath*p; rhs];
    col = vecnorm(A); col(col==0)=1;
    A_s = A ./ col;
    dp_s = (A_s'*A_s + lambda*eye(K)) \ (A_s'*r);
    dp = dp_s ./ col';
    p = p + dp;
    if norm(dp) < 1e-6, break; end
end

dg = p(1:Kint); Gext = p(Kint+1:end);
A = [Apath; H];
col = vecnorm(A); col(col==0)=1;
A_s = A ./ col;
info.cond = cond(A_s);
info.svd = svd(A_s);
info.iterations = it;
info.Cov_s = inv(A_s'*A_s + lambda*eye(K));
end
