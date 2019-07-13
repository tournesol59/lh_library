void demo_mexfun(double * out, const double * in, const int size) {
  
   int i;
   for (i=0; i<size;i++) {
      out[i]=2*in[i];
   }
}

