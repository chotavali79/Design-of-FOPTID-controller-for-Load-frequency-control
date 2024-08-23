function [z] = levy_FOPID(pop,m,omega)  
    num = gamma(1+omega)*sin(pi*omega/2); % used for Numerator 
    den = gamma((1+omega)/2)*omega*2^((omega-1)/2); % used for Denominator
    sigma_u = (num/den)^(1/omega);% Standard deviation
    u = random('Normal',0,sigma_u,pop,m); 
    v = random('Normal',0,1,pop,m);
    z =u./(abs(v).^(1/omega));
  end