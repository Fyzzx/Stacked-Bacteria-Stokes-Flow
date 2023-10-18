function yyy = orthoNorm(xxx)
    xxxmag = sqrt(xxx(:,1).^2 + xxx(:,2).^2 + xxx(:,3).^2);
    yyy(:,1) = xxx(:,1)./xxxmag; 
    yyy(:,2) = xxx(:,2)./xxxmag;
    yyy(:,3) = xxx(:,3)./xxxmag;
end