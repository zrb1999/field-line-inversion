function [dg_err,B_err,L_err]=ensemble_metrics(obs,dg_true,dg,f,R_E)
dg_err=norm(dg-dg_true)/norm(dg_true);
Bs=[]; rr=[];
for j=1:numel(obs)
  o=obs(j);
  for i=1:numel(o.r_true)
    [~,~,Bt]=field_from_flux(o.fun,o.r_true(i),o.theta(i));
    [~,~,Be]=field_from_flux(f,o.r_true(i),o.theta(i));
    Bs=[Bs;Be-Bt]; %#ok<AGROW>
  end
  C=f(o.r_obs(1),o.theta(1));
  re=zeros(numel(o.r_true),1);
  re(1)=o.r_obs(1);
  for i=2:numel(o.r_true)
      % Trace the reconstructed level set from the anchor / previous point,
      % not from the noisy observed radius at each point.
      re(i)=trace_level_set_fun(o.theta(i),C,f,re(i-1),[0.9,3.4]);
  end
  rr=[rr;re-o.r_true]; %#ok<AGROW>
end
B_err=sqrt(mean(Bs.^2));
L_err=sqrt(mean(rr.^2))*R_E;
end
