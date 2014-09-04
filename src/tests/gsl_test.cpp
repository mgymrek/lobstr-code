#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_cblas.h>
#include <err.h>

int
main (void)
{
  double x = 5.0;
  double y = gsl_sf_bessel_J0 (x);
  printf ("J0(%g) = %.18e\n", x, y);

  /* Expected result, as in
     https://www.gnu.org/software/gsl/manual/html_node/An-Example-Program.html#An-Example-Program

     -1.775967713143382920e-01
  */
  if ( y < -0.177597 || y>-0.177596 )
    errx(1,"gsl calculation test failed (y=%g)", y);
  return 0;
}
