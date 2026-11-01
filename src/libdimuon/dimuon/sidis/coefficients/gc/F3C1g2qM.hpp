const double c7t8 = p.c.nlo.c7t8;
const double c8t8 = p.c.nlo.c8t8;
const double c9t8 = p.c.nlo.c9t8;
return ((c7t8*math::power<-1>((1.+(-1.*p.z)))*p.ff*p.pdf_hat)+(c8t8*p.ff_hat*p.pdf_hat)+(c9t8*((math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+(math::power<-1>((1.+(-1.*p.z)))*p.ff*std::log((1.+(-1.*p.z)))))*p.pdf_hat));