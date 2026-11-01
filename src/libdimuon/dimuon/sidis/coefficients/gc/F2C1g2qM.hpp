const double c7t7 = p.c.nlo.c7t7;
const double c8t7 = p.c.nlo.c8t7;
const double c9t7 = p.c.nlo.c9t7;
return ((c7t7*math::power<-1>((1.+(-1.*p.z)))*p.ff*p.pdf_hat)+(c8t7*p.ff_hat*p.pdf_hat)+(c9t7*((math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+(math::power<-1>((1.+(-1.*p.z)))*p.ff*std::log((1.+(-1.*p.z)))))*p.pdf_hat));