function [dg, Gext, Cov, info] = invert_flux_perturbation_full(obs, g_base, Ni, Ne, lambda)
% Invert perturbations relative to a known internal baseline g_base.
% Full internal perturbation vector: dg_1..dg_Ni (including the dipole).
% Unknowns: [dg_1..dg_Ni, G_1..G_Ne].
%
% This is the variant used for the IGRF secular-variation experiments, where
% the dipole also changes between epochs.  The original
% invert_flux_perturbation.m fixes g10 and solves only dg_2..dg_Ni.

Kint = Ni; Kext = Ne; K = Kint + Kext;
A = []; b = [];

for j = 1:numel(obs)
    r_obs = obs(j).r_obs;
    theta_obs = obs(j).theta;
    M = numel(r_obs) - 1;

    r_a = r_obs(1); th_a = theta_obs(1);
    Xi0_a = 0;
    for n = 1:Ni
        [Phi,~,~] = flux_basis(r_a, th_a, n);
        Xi0_a = Xi0_a + g_base(n)*Phi;
    end

    Aj = zeros(M,K);
    bj = zeros(M,1);
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
    A = [A; Aj]; %#ok<AGROW>
    b = [b; bj]; %#ok<AGROW>
end

col = vecnorm(A); col(col==0)=1;
A_s = A ./ col;
H = A_s'*A_s + lambda*eye(K);
p_s = H \ (A_s'*b);
p = p_s ./ col';

dg = p(1:Kint);
Gext = p(Kint+1:end);
Cov_s = inv(H);
Cov = diag(1./col)*Cov_s*diag(1./col);
info.cond = cond(A_s);
info.svd = svd(A_s);
info.colnorm = col;
end
