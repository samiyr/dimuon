const double c7t24 = p.c.nnlo.c7t24;
const double c8t24 = p.c.nnlo.c8t24;
const double c9t24 = p.c.nnlo.c9t24;
const double c10t24 = p.c.nnlo.c10t24;
const double c11t24 = p.c.nnlo.c11t24;
return ((c7t24*math::power<-1>((1.+(-1.*p.z)))*p.ff*p.pdf_hat)+(c8t24*p.ff_hat*p.pdf_hat)+(c9t24*((math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+(math::power<-1>((1.+(-1.*p.z)))*p.ff*std::log((1.+(-1.*p.z)))))*p.pdf_hat)+(c10t24*((logv1*math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+((1./2.)*math::power<-1>((1.+(-1.*p.z)))*p.ff*math::power<2>(std::log((1.+(-1.*p.z))))))*p.pdf_hat)+(c11t24*((math::power<2>(logv1)*math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+((1./3.)*math::power<-1>((1.+(-1.*p.z)))*p.ff*math::power<3>(std::log((1.+(-1.*p.z))))))*p.pdf_hat));