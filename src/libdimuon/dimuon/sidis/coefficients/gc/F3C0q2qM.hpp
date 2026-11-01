const double c1t2 = p.c.lo.c1t2;
return (c1t2*math::power<-1>((1.+(-1.*p.pdf_momentum_fraction)))*math::power<-1>((1.+(-1.*p.z)))*p.ff*p.pdf);