const double c7t6 = p.c.nlo.c7t6;
const double c8t6 = p.c.nlo.c8t6;
const double c9t6 = p.c.nlo.c9t6;
return ((c7t6*math::power<-1>((1.+(-1.*p.z)))*p.ff*p.pdf_hat)+(c8t6*p.ff_hat*p.pdf_hat)+(c9t6*((math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+(math::power<-1>((1.+(-1.*p.z)))*p.ff*std::log((1.+(-1.*p.z)))))*p.pdf_hat));