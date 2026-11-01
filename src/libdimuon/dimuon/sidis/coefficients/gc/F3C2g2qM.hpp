const double c7t25 = p.c.nnlo.c7t25;
const double c8t25 = p.c.nnlo.c8t25;
const double c9t25 = p.c.nnlo.c9t25;
const double c10t25 = p.c.nnlo.c10t25;
const double c11t25 = p.c.nnlo.c11t25;
return ((c7t25*math::power<-1>((1.+(-1.*p.z)))*p.ff*p.pdf_hat)+(c8t25*p.ff_hat*p.pdf_hat)+(c9t25*((math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+(math::power<-1>((1.+(-1.*p.z)))*p.ff*std::log((1.+(-1.*p.z)))))*p.pdf_hat)+(c10t25*((logv1*math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+((1./2.)*math::power<-1>((1.+(-1.*p.z)))*p.ff*math::power<2>(std::log((1.+(-1.*p.z))))))*p.pdf_hat)+(c11t25*((math::power<2>(logv1)*math::power<-1>((1.+(-1.*p.zeta)))*((-1.*p.ff)+p.ff_hat))+((1./3.)*math::power<-1>((1.+(-1.*p.z)))*p.ff*math::power<3>(std::log((1.+(-1.*p.z))))))*p.pdf_hat));