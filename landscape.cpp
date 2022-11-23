//The program of TDA
// by Xu Fei 2021.4.24
#include <fstream>
#include <iostream>
#include <typeinfo> 
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <cstring>
#include <vector>

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

/* Period parameters */  
#define NNNN 624
#define MMMM 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UMASK 0x80000000UL /* most significant w-r bits */
#define LMASK 0x7fffffffUL /* least significant r bits */
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))

static unsigned long state[NNNN]; /* the array for the state vector  */
static int left = 1;
static int initf = 0;
static unsigned long *next;

/* initializes state[NNNN] with a seed */
void init_genrand(unsigned long s)
{
    int j;
    state[0]= s & 0xffffffffUL;
    for (j=1; j<NNNN; j++) {
        state[j] = (1812433253UL * (state[j-1] ^ (state[j-1] >> 30)) + j); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array state[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        state[j] &= 0xffffffffUL;  /* for >32 bit machines */
    }
    left = 1; initf = 1;
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array(unsigned long init_key[], unsigned long key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (NNNN>key_length ? NNNN : key_length);
    for (; k; k--) {
        state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1664525UL)) + init_key[j] + j; /* non linear */
        state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=NNNN) { state[0] = state[NNNN-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=NNNN-1; k; k--) {
        state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=NNNN) { state[0] = state[NNNN-1]; i=1; }
    }

    state[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
    left = 1; initf = 1;
}

static void next_state(void)
{
    unsigned long *p=state;
    int j;

    /* if init_genrand() has not been called, */
    /* a default initial seed is used         */
    if (initf==0) init_genrand(5489UL);

    left = NNNN;
    next = state;
    
    for (j=NNNN-MMMM+1; --j; p++) 
        *p = p[MMMM] ^ TWIST(p[0], p[1]);

    for (j=MMMM; --j; p++) 
        *p = p[MMMM-NNNN] ^ TWIST(p[0], p[1]);

    *p = p[MMMM-NNNN] ^ TWIST(p[0], state[0]);
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;

    if (--left == 0) next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    unsigned long y;

    if (--left == 0) next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (long)(y>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    unsigned long y;

    if (--left == 0) next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (double)y * (1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    unsigned long y;

    if (--left == 0) next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (double)y * (1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    unsigned long y;

    if (--left == 0) next_state();
    y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return ((double)y + 0.5) * (1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

#define NT (int)1.44e6
#define dt (double)0.001
#define NP (double)10000
#define pi (double)3.1415926

int i, j, it, ii, jj, kk;

double TRADD, RIPK1, Casp8, RIPK3, Acasp8, sumC8;
double temp_TRADD, temp_RIPK1, temp_Casp8, temp_RIPK3, temp_Acasp8, temp_sumC8;

double V1, V2, V3, V4, V5, V6;

double const k0=1.7;
double const k1=9.5;
double const k2=0.35;
double const k3=6.7;
double const k4=0.17;
double const k5=1.4;
double const k6=0.3; //0.3;
double const k7=0.0;
double const k8=0.2;
double const k9=2.1;
double const k10=1.0;
double const k11=8.5;
double const J0=0.012;
double const J1=0.12;
double const J2=2.3;
double const J3=1.2;
double const J4=1.47;
double const J5=0.008;
double const J6=10.4;
double const J7=44.2;
double const J8=0.036;
double const J9=0.16;
double const J10=11.4;
double const J11=0.0015;
double const d0=0.03;
double const d1=0.1;
double const d2=0.35;
double const d3=0.14;
double const n0=3.0;
double const n1=2.0;
double const n2=1.0;
double const n3=3.0;
double const n4=4.0;
double const n5=4.0;
double const n6=2.0;
double const n7=2.0;
double const n8=4.0;
double const n9=4.0;
double const n10=4.0;
double const n11=4.0;
double const TNF=0.5;

double const k12=3.6;
double const J12=1.25;
double const n12=4.0;
double const d4=0.35;

double const g=0.115;

int size_space = 200;

double rand_numa1, rand_numb1, rand_numa2, rand_numb2, rand_numa3, rand_numb3, rand_numa4, rand_numb4, rand_numa5, rand_numb5, rand_numa6, rand_numb6;

double t, n_sth;

int main()
{
	std::vector< std::vector<double> > paras(NP, std::vector<double>(5, 0.0));
	std::ifstream infile;
	infile.open("LHS_of_paras.dat");
	for(i=0;i<NP;i++)
	{
		for(j=0;j<5;j++)
		{
			infile >> paras[i][j];
		}
	}
	infile.close();
	
	std::vector< std::vector<double> > phase_space(size_space, std::vector<double>(size_space, 0.0));
	
	std::ofstream output_file("trajector_density.dat");
	
	n_sth=5.0e-3;
	
	srand((unsigned)time(NULL));
	unsigned long init[4]={rand(), rand(), rand(), rand()},  length=4;
	//unsigned long init[4]={0x283, 0x175, 0x615, 0x954},length=4;
	init_by_array(init, length);
	
	for(i=0;i<NP;i++)
	{
		TRADD=paras[i][0];
		RIPK1=g*paras[i][1];
		RIPK3=paras[i][2];
		sumC8=paras[i][3];
		Casp8=sumC8*1e-5;
		Acasp8=sumC8*(1-1e-5);
		
		printf("===i=%d====\n",i);
		//printf("===i=%d==TRADD=%f==RIP1=%f==Casp8=%f==RIP3=%f==\n",i,TRADD,RIPK1,Casp8,RIPK3);
		
		int S_TRADD = 0;
        int S_RIPK1 = 0;
        int S_Casp8 = 0;
        int S_RIPK3 = 0;
        int S_Acasp8 = 0;
        int S_sumC8 = 0;
        
		for(it=1; it<=NT; it++)
		{
			rand_numa1 = genrand_real1(); rand_numb1 = genrand_real1();
			rand_numa2 = genrand_real1(); rand_numb2 = genrand_real1();
			rand_numa3 = genrand_real1(); rand_numb3 = genrand_real1();
			rand_numa4 = genrand_real1(); rand_numb4 = genrand_real1();
			rand_numa5 = genrand_real1(); rand_numb5 = genrand_real1();
			rand_numa6 = genrand_real1(); rand_numb6 = genrand_real1();
			
			temp_TRADD=TRADD; temp_RIPK1=RIPK1; temp_Casp8=Casp8; temp_RIPK3=RIPK3, temp_Acasp8=Acasp8, temp_sumC8=sumC8;
			
			t=it*dt;
			
			V1=k0*(1-TRADD)*pow(TNF,n0)/(pow(TNF,n0)+pow(J0,n0)) - k1*TRADD*pow(RIPK1,n1)/(pow(RIPK1,n1)+pow(J1,n1)) - d0*TRADD;
			V2=k2*(g-RIPK1)*pow(TNF,n2)/(pow(TNF,n2)+pow(J2,n2)) + k3*(g-RIPK1)*pow(RIPK3,n3)/(pow(RIPK3,n3)+pow(J3,n3)) - k4*RIPK1*pow(TRADD,n4)/(pow(TRADD,n4)+pow(J4,n4)) - k5*RIPK1*pow(Casp8,n5)/(pow(Casp8,n5)+pow(J5,n5)) - d1*RIPK1;
			V3=k6*(1-Acasp8-Casp8)*pow(RIPK1,n6)/(pow(RIPK1,n6)+pow(J6,n6)) + k7*(1-Casp8)*pow(TRADD,n7)/(pow(TRADD,n7)+pow(J7,n7)) - k8*Casp8*pow(RIPK3,n8)/(pow(RIPK3,n8)+pow(J8,n8)) - d2*Casp8;
			V4=k9*(1-RIPK3)*pow(RIPK1,n9)/(pow(RIPK1,n9)+pow(J9,n9)) + k10*(1-RIPK3)*pow(RIPK3,n10)/(pow(RIPK3,n10)+pow(J10,n10)) - k11*RIPK3*pow(Casp8,n11)/(pow(Casp8,n11)+pow(J11,n11)) - d3*RIPK3;
			V5=k12*(1-Acasp8-Casp8)*pow(TRADD,n12)/(pow(TRADD,n12)+pow(J12,n12))-d4*Acasp8;
			V6=V3+V5;
			
			TRADD=temp_TRADD + dt*V1 + temp_TRADD*sqrt(2*n_sth*dt)*(0.0 + (sqrt(-2.0*log(rand_numa1)) * cos(2.0*pi*rand_numb1))*1.0);
			
			RIPK1=temp_RIPK1 + dt*V2 + temp_RIPK1*sqrt(2*n_sth*dt)*(0.0 + (sqrt(-2.0*log(rand_numa2)) * cos(2.0*pi*rand_numb2))*1.0);
			
			Casp8=temp_Casp8 + dt*V3 + temp_Casp8*sqrt(2*n_sth*dt)*(0.0 + (sqrt(-2.0*log(rand_numa3)) * cos(2.0*pi*rand_numb3))*1.0);
			
			RIPK3=temp_RIPK3 + dt*V4 + temp_RIPK3*sqrt(2*n_sth*dt)*(0.0 + (sqrt(-2.0*log(rand_numa4)) * cos(2.0*pi*rand_numb4))*1.0);
			
			Acasp8=temp_Acasp8 + dt*V5 + temp_Acasp8*sqrt(2*n_sth*dt)*(0.0 + (sqrt(-2.0*log(rand_numa5)) * cos(2.0*pi*rand_numb5))*1.0);
			
			sumC8=temp_sumC8 + dt*V6 + temp_Casp8*sqrt(2*n_sth*dt)*(0.0 + (sqrt(-2.0*log(rand_numa3)) * cos(2.0*pi*rand_numb3))*1.0) + sqrt(2*temp_Acasp8*n_sth*dt)*(0.0 + (sqrt(-2.0*log(rand_numa5)) * cos(2.0*pi*rand_numb5))*1.0);
					
			if(TRADD<0.0)
			{
				TRADD = 0.0;
			}
			if(TRADD>1.0)
			{
				TRADD = 1.0;
			}
			if(RIPK1<0.0)
			{
				RIPK1 = 0.0;
			}
			if(RIPK1>g)
			{
				RIPK1 = g;
			} 
			if(Casp8<0.0)
			{
				Casp8 = 0.0;
			}
			if(Casp8>1.0)
			{
				Casp8 = 1.0;
			} 
			if(RIPK3<0.0)
			{
				RIPK3 = 0.0;
			}
			if(RIPK3>1.0)
			{
				RIPK3 = 1.0;
			}
			if(Acasp8<0.0)
			{
				Acasp8 = 0.0;
			}
			if(Acasp8>1.0)
			{
				Acasp8 = 1.0;
			}
			if(sumC8<0.0)
			{
				sumC8 = 0.0;
			}
			if(sumC8>1.0)
			{
				sumC8 = 1.0;
			}
			
			if(it%10==0)
			{
				if(it==NT)
				{
					printf("====RIPK3=%lf====sumcasp8=%lf====\n",RIPK3,sumC8);
				}
				
				std::vector<int> space_position(6, 0);
				
				for(kk=0;kk<size_space;kk++)
				{
					if(TRADD < (kk+1)*1.0/size_space)
					{
						space_position[0] = kk;
                        break;
					}
				}
				
				for(kk=0;kk<size_space;kk++)
				{
					if(RIPK1 < (kk+1)*1.0/size_space)
					{
						space_position[1] = kk;
                        break;
					}
				}
				
				for(kk=0;kk<size_space;kk++)
				{
					if(Casp8 < (kk+1)*1.0/size_space)
					{
						space_position[2] = kk;
                        break;
					}
				}
				
				for(kk=0;kk<size_space;kk++)
				{
					if(RIPK3 < (kk+1)*1.0/size_space)
					{
						space_position[3] = kk;
                        break;
					}
				}
				
				for(kk=0;kk<size_space;kk++)
				{
					if(Acasp8 < (kk+1)*1.0/size_space)
					{
						space_position[4] = kk;
                        break;
					}
				}
				
				for(kk=0;kk<size_space;kk++)
				{
					if(sumC8 < (kk+1)*1.0/size_space)
					{
						space_position[5] = kk;
                        break;
					}
				}
				
				if(space_position[3] != S_RIPK3 || space_position[5] != S_sumC8)
				{
	                S_RIPK3 = space_position[3];
					S_sumC8 = space_position[5];
	                //printf("====S_RIPK3=%d====S_sumC8=%d====\n",S_RIPK3,S_sumC8);
	                phase_space[S_RIPK3][S_sumC8] += 1.0;
	            }
	            
			}
						
		}
		
	}
	double sum_PS=0.0;
	double U;
	for(i=0;i<size_space;i++)
	{
		for(j=0;j<size_space;j++)
		{
			sum_PS += phase_space[i][j];
		}
	}
	for(i=0;i<size_space;i++)
	{
		for(j=0;j<size_space;j++)
		{
			if(phase_space[i][j]==0)
			{
				U = 30.0;
			}
			else
			{
				U = -log(phase_space[i][j]/sum_PS);
				if(U>30.0)
				{
					U = 30.0;
				}
			}
			char out_char[100];
			sprintf(out_char, "%f\t%f\t%f\n",((i+1.0)*1.0/size_space),((j+1.0)*1.0/size_space),U);
			output_file << out_char;
		}
		//output_file << std::endl;
	}
	
	output_file.close();
	return 0;
}
