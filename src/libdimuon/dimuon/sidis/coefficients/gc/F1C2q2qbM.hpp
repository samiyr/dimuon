const double c2t8 = p.c.nnlo.c2t8;
const double c7t8 = p.c.nnlo.c7t8;
const double c8t8 = p.c.nnlo.c8t8;
const double c9t8 = p.c.nnlo.c9t8;
const double c14t8 = p.c.nnlo.c14t8;
return ((c2t8*math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*p.ff_hat*p.pdf)+(c7t8*math::power<-1>((1.+(-1.*p.z)))*p.ff*p.pdf_hat)+(c8t8*p.ff_hat*p.pdf_hat)+(c9t8*((math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+(math::power<-1>((1.+(-1.*p.z)))*p.ff*std::log((1.+(-1.*p.z)))))*p.pdf_hat)+(c14t8*p.ff_hat*((math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*std::log((1.+(-1.*p.pdf_momentum_fraction)))*p.pdf)+(math::power<-1>((1.+(-1.*p.xi)))*((-1.*p.pdf)+p.pdf_hat)))));