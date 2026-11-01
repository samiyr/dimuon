const double c2t10 = p.c.nnlo.c2t10;
const double c7t10 = p.c.nnlo.c7t10;
const double c8t10 = p.c.nnlo.c8t10;
const double c9t10 = p.c.nnlo.c9t10;
const double c14t10 = p.c.nnlo.c14t10;
return ((c2t10*math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*p.ff_hat*p.pdf)+(c7t10*math::power<-1>((1.+(-1.*p.z)))*p.ff*p.pdf_hat)+(c8t10*p.ff_hat*p.pdf_hat)+(c9t10*((math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+(math::power<-1>((1.+(-1.*p.z)))*p.ff*std::log((1.+(-1.*p.z)))))*p.pdf_hat)+(c14t10*p.ff_hat*((math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*std::log((1.+(-1.*p.pdf_momentum_fraction)))*p.pdf)+(math::power<-1>((1.+(-1.*p.xi)))*((-1.*p.pdf)+p.pdf_hat)))));