const double c7t23 = p.c.nnlo.c7t23;
const double c8t23 = p.c.nnlo.c8t23;
const double c9t23 = p.c.nnlo.c9t23;
const double c10t23 = p.c.nnlo.c10t23;
const double c11t23 = p.c.nnlo.c11t23;
return ((c7t23*math::power<-1>((1.+(-1.*p.z)))*p.ff*p.pdf_hat)+(c8t23*p.ff_hat*p.pdf_hat)+(c9t23*((math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+(math::power<-1>((1.+(-1.*p.z)))*p.ff*std::log((1.+(-1.*p.z)))))*p.pdf_hat)+(c10t23*((logv1*math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+((1./2.)*math::power<-1>((1.+(-1.*p.z)))*p.ff*math::power<2>(std::log((1.+(-1.*p.z))))))*p.pdf_hat)+(c11t23*((math::power<2>(logv1)*math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+((1./3.)*math::power<-1>((1.+(-1.*p.z)))*p.ff*math::power<3>(std::log((1.+(-1.*p.z))))))*p.pdf_hat));