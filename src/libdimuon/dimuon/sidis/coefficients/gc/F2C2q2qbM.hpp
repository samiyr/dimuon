const double c2t9 = p.c.nnlo.c2t9;
const double c7t9 = p.c.nnlo.c7t9;
const double c8t9 = p.c.nnlo.c8t9;
const double c9t9 = p.c.nnlo.c9t9;
const double c14t9 = p.c.nnlo.c14t9;
return ((c2t9*math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*p.ff_hat*p.pdf)+(c7t9*math::power<-1>((1.+(-1.*p.z)))*p.ff*p.pdf_hat)+(c8t9*p.ff_hat*p.pdf_hat)+(c9t9*((math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+(math::power<-1>((1.+(-1.*p.z)))*p.ff*std::log((1.+(-1.*p.z)))))*p.pdf_hat)+(c14t9*p.ff_hat*((math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*std::log((1.+(-1.*p.pdf_momentum_fraction)))*p.pdf)+(math::power<-1>((1.+(-1.*p.xi)))*((-1.*p.pdf)+p.pdf_hat)))));