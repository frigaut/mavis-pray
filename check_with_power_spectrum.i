func filt_power_spectrum(D,slo)
{
  // D = 8.;
  // slo = -2.5;
  sfv = span(0.01,100,1000); // spatial frequencies vector in cycle-1
  psp = sfv^slo;
  fil = (1.-2*bessj1(pi*D*sfv)/(pi*D*sfv));
  fil = sfv>(1./D);
  plot,psp,sfv;
  plg,psp*fil,sfv,color="red";
  logxy,1,1;
  pow = sum(psp*sfv);
  powf = sum(psp*sfv*fil);
  write,format="Residual = %.2f%%\n",100.*powf/pow;
  return powf/pow;
}
