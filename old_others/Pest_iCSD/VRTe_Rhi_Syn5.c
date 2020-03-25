#include <stdio.h>
#include <gsl/gsl_vector.h>
/* vf vector first, vs vector second
 * ff file first, fs file second
 * sf size first, ss size second*/

int
main (void)
{
  int i;
  int m;
  int n;
  int vf_index;
  int sf = 1020;  /* Totsl VRTe simulation */
  int ss = 5; /* num weight, also = num VRTe */
  int so = 204; /* VRTe linear combinatin according to weight */

  gsl_vector * vf = gsl_vector_alloc (sf);
  gsl_vector * vs = gsl_vector_alloc (ss);
  gsl_vector * out = gsl_vector_alloc (so);
  gsl_vector_set_zero (out);
 
  {
  FILE * ff = fopen ("VRTe.sim","r");
  gsl_vector_fscanf (ff, vf);
  fclose (ff);
  }
  {
  FILE * fs = fopen ("Weight.wgt","r");
  gsl_vector_fscanf (fs, vs);
  fclose (fs);
  }

  /* sum */
  for (m = 0; m < ss; m++) /* weights */
  {
    for (n = 0; n < so; n++) /* VRTeSimulations */
    {
      vf_index = n + m * so;
      /*printf("vf_index %g\n",vf_index);*/ 
      gsl_vector_set(out, n, gsl_vector_get(out, n) + gsl_vector_get(vf, vf_index) * gsl_vector_get(vs, m));
    }
  }  
 
  FILE * fo = fopen("PESTinput.txt","w");
  fputs("# data\n",fo);
  gsl_vector_fprintf(fo,out,"%.10f");
  fclose(fo);

  gsl_vector_free (vf);
  gsl_vector_free (vs);
  gsl_vector_free (out);

  return 0;
}
