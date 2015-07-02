#include <stdio.h>
#include <math.h>
#include <fenv.h>

int main(void)
{
    double ans = 1., tmp = 1;
    double eps = 1;
    double x = -19.5;
    long count = 0;

    printf("rounding mode = %d\n", fegetround());
    fesetround(FE_UPWARD);
    printf("rounding mode = %d\n", fegetround());
   
    while (eps > 1e-15) {
        count++;
        tmp *= x/(count);
        eps = fabs(tmp);
        ans += tmp;
//        printf("count = %d, eps = %.8le\n", count, eps);
    }
    
    printf("math.h exp(%.6lf) = %.9le\n", x, exp(x));
    printf("ans = %.16le\n", ans);

    return 0;
}

