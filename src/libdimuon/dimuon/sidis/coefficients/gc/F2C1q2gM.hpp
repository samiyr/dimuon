const double c2t4 = p.c.nlo.c2t4;
const double c8t4 = p.c.nlo.c8t4;
const double c14t4 = p.c.nlo.c14t4;
return ((c2t4*math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*p.ff_hat*p.pdf)+(c8t4*p.ff_hat*p.pdf_hat)+(c14t4*p.ff_hat*((math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*std::log((1.+(-1.*p.pdf_momentum_fraction)))*p.pdf)+(math::power<-1>((1.+(-1.*p.xi)))*((-1.*p.pdf)+p.pdf_hat)))));