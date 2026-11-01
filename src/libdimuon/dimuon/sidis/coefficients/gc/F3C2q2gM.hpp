const double c2t22 = p.c.nnlo.c2t22;
const double c8t22 = p.c.nnlo.c8t22;
const double c14t22 = p.c.nnlo.c14t22;
const double c20t22 = p.c.nnlo.c20t22;
const double c26t22 = p.c.nnlo.c26t22;
return ((c2t22*math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*p.ff_hat*p.pdf)+(c8t22*p.ff_hat*p.pdf_hat)+(c14t22*p.ff_hat*((math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*std::log((1.+(-1.*p.pdf_momentum_fraction)))*p.pdf)+(math::power<-1>((1.+(-1.*p.xi)))*((-1.*p.pdf)+p.pdf_hat))))+(c20t22*p.ff_hat*(((1./2.)*math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*math::power<2>(std::log((1.+(-1.*p.pdf_momentum_fraction))))*p.pdf)+(logv3*math::power<-1>((1.+(-1.*p.xi)))*((-1.*p.pdf)+p.pdf_hat))))+(c26t22*p.ff_hat*(((1./3.)*math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*math::power<3>(std::log((1.+(-1.*p.pdf_momentum_fraction))))*p.pdf)+(math::power<2>(logv3)*math::power<-1>((1.+(-1.*p.xi)))*((-1.*p.pdf)+p.pdf_hat)))));