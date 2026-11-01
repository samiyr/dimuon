const double c2t14 = p.c.nnlo.c2t14;
const double c8t14 = p.c.nnlo.c8t14;
const double c14t14 = p.c.nnlo.c14t14;
const double c20t14 = p.c.nnlo.c20t14;
return ((c2t14*math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*p.ff_hat*p.pdf)+(c8t14*p.ff_hat*p.pdf_hat)+(c14t14*p.ff_hat*((math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*std::log((1.+(-1.*p.pdf_momentum_fraction)))*p.pdf)+(math::power<-1>((1.+(-1.*p.xi)))*((-1.*p.pdf)+p.pdf_hat))))+(c20t14*p.ff_hat*(((1./2.)*math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*math::power<2>(std::log((1.+(-1.*p.pdf_momentum_fraction))))*p.pdf)+(logv3*math::power<-1>((1.+(-1.*p.xi)))*((-1.*p.pdf)+p.pdf_hat)))));