const double c2t3 = p.c.nlo.c2t3;
const double c8t3 = p.c.nlo.c8t3;
const double c14t3 = p.c.nlo.c14t3;
return ((c2t3*math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*p.ff_hat*p.pdf)+(c8t3*p.ff_hat*p.pdf_hat)+(c14t3*p.ff_hat*((math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*std::log((1.+(-1.*p.pdf_momentum_fraction)))*p.pdf)+(math::power<-1>((1.+(-1.*p.xi)))*((-1.*p.pdf)+p.pdf_hat)))));