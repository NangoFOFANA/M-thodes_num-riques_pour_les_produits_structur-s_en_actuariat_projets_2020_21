%%% fonction g(x)
function G =fonc(x)
       mu_i=0.04;
       eta=2;
       G=mu_i+eta*max(-x,0);
