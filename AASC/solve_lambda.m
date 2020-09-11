function F = solve_lambda(alpha,gamma,ini)
    %fsolve
    k=size(gamma,2);
    a=ones(1,k);
    f=@(x)(abs((sum(a./((gamma-x).*alpha))^2)/(sum(a./(((gamma-x).^2).*alpha))))-1);
    %initial x0
    if ini==0
        x0=-1;
    elseif ini==k+1
        x0=1;
    else
        x0=round(gamma(1,ini)*100)/100;
    end
    F=fsolve(f,x0, optimoptions('fsolve','Display','none'));
    