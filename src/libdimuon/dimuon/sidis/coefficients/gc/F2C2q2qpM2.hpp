const double c7t18 = p.c.nnlo.c7t18;
const double c8t18 = p.c.nnlo.c8t18;
const double c9t18 = p.c.nnlo.c9t18;
const double c10t18 = p.c.nnlo.c10t18;
return ((c7t18*math::power<-1>((1.+(-1.*p.z)))*p.ff*p.pdf_hat)+(c8t18*p.ff_hat*p.pdf_hat)+(c9t18*((math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+(math::power<-1>((1.+(-1.*p.z)))*p.ff*std::log((1.+(-1.*p.z)))))*p.pdf_hat)+(c10t18*((logv1*math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+((1./2.)*math::power<-1>((1.+(-1.*p.z)))*p.ff*math::power<2>(std::log((1.+(-1.*p.z))))))*p.pdf_hat));