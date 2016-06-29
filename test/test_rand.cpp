


#include "../lib/aos.h"

rand_t rstat;/*a struct for each sequence*/
rand_t rstat2;


static void test_number(){
    rand_t a,b;
    seed_rand(&a,1);
    info("First number is %g\n",randn(&a));
    seed_rand(&a,1);
    unsigned int i=lrand(&a);
    info("First number is %u\n",i);
    i=lrand(&a);
    info("First number is %u\n",i);
    seed_rand(&b,i);
    info("First number is %g\n",randn(&b));
}
static void test_speed(){
    double *p;
    int N,i;
    clock_t ck;
    int seed=2;
    N=1000000;
    p=mymalloc(N,double);
    memset(p, 0, sizeof(double)*N);
    seed_rand(&rstat, seed);
    srand(seed);
    ck=clock();
    seed_rand(&rstat, seed);
    ck=clock();
    for(i=0; i<N; i++){
	p[i]=randu(&rstat);
    }
    printf("randu elapsed %f seconds\n", (double)(clock()-ck)/CLOCKS_PER_SEC);
    writedbl(p, N, 1,"randu.bin");
    srand(1);
    ck=clock();
    for(i=0; i<N; i++){
	p[i]=rand();
    }
    printf("Rand elapsed %f seconds\n", (double)(clock()-ck)/CLOCKS_PER_SEC);
    /*ck=clock();
    for(i=0; i<N; i++){
	p[i]=slow_randn(&rstat);
    }
    printf("Slow Randn elapsed %f seconds\n", (double)(clock()-ck)/CLOCKS_PER_SEC);
    writedbl("slow_randn.bin", p, N, 1);
    */
    ck=clock();
    for(i=0; i<N; i++){
	p[i]=randn(&rstat)*10;
    }
    printf("Fast Randn elapsed %f seconds\n", (double)(clock()-ck)/CLOCKS_PER_SEC);
    writedbl(p, N, 1, "randn.bin");

    seed_rand(&rstat2,seed);
    ck=clock();
    for(i=0; i<N; i++){
	p[i]=(double)randp(&rstat2,30);
    }
    printf("Randp elapsed %f seconds\n", (double)(clock()-ck)/CLOCKS_PER_SEC);
    writedbl(p, N, 1,"randp.bin");
    free(p);
}
int main(){
    test_number();
    test_speed();
    return 0;
}
