const double c7t19 = p.c.nnlo.c7t19;
const double c8t19 = p.c.nnlo.c8t19;
const double c9t19 = p.c.nnlo.c9t19;
const double c10t19 = p.c.nnlo.c10t19;
return ((c7t19*math::power<-1>((1.+(-1.*p.z)))*p.ff*p.pdf_hat)+(c8t19*p.ff_hat*p.pdf_hat)+(c9t19*((math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+(math::power<-1>((1.+(-1.*p.z)))*p.ff*std::log((1.+(-1.*p.z)))))*p.pdf_hat)+(c10t19*((logv1*math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+((1./2.)*math::power<-1>((1.+(-1.*p.z)))*p.ff*math::power<2>(std::log((1.+(-1.*p.z))))))*p.pdf_hat));