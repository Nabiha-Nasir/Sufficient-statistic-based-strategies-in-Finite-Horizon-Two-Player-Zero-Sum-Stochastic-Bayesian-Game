function [p_present] = update_belief_li_pplus(k,p, X_star,a,b,P)
    X_star(isnan(X_star))=0;
    p_temp=p.*X_star(a,:);
    xbar=sum(p_temp);
    if xbar<10^(-7)
        p_present=zeros(1,k);
    else
        p_present=p_temp*P{a,b}/xbar;
    end
   