//The program of soucan
// by Xu Fei 2021.10.18
#include <fstream>
#include <iostream>
#include <typeinfo> 
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <cstring>
#include <vector>

#include"mpi.h"

#define NT (int)7.2e6
#define dt (double)0.0001
#define NP (int)10000

int i, j, it, ii, jj, kk;

double TNF, TRADD, RIPK1, Casp8, RIPK3, aCasp8;
double temp_TRADD, temp_RIPK1, temp_Casp8, temp_RIPK3, temp_aCasp8; 

double g, t;

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);//MPI库的初始化
	int numprocs, myid;//定义进程总数及进程标识
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);//获得进程数
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);//获得当前进程标识号0,1,2,3,....,numprocs - 1
	
        //start = MPI_Wtime();//获得时间
//============================================================================================================= 
	
	float paras[NP][29]={0.0};
	std::ifstream infile;
	infile.open("LHS_of_paras.dat");
	for(i=0;i<NP;i++)
	{
		for(j=0;j<29;j++)
		{
			infile >> paras[i][j];
		}
	}
	infile.close();
	
	float Hills[NP][12]={0.0};
	std::ifstream Hillfile;
	Hillfile.open("Hill.dat");
	for(i=0;i<NP;i++)
	{
		for(j=0;j<12;j++)
		{
			Hillfile >> Hills[i][j];
		}
	}
	Hillfile.close();
	
	std::ofstream rise_file("rise.dat",std::ofstream::out| std::ofstream::app); //| std::ofstream::ate
	std::ofstream Hvalue_file("h_value.dat",std::ofstream::out| std::ofstream::app);
	
	std::vector<double> k(12,0.0);
	std::vector<double> J(12,0.0);	
	std::vector<double> n(12,0.0);	
	std::vector<double> d(4,0.0);	
	
	std::vector<double> RIPK3_2h(50,0.0);
	
	int maxindex;
	
	double max_RIPK3, last_RIPK3;
	
	for(i = myid; i < NP; i += numprocs)
	{
		for(jj=0;jj<12;jj++)
		{
			k[jj] = paras[i][jj];
			J[jj] = paras[i][jj+12];
		}
		
		for(jj=0;jj<4;jj++)
		{
			d[jj] = paras[i][jj+24];
		}
		
		TNF=paras[i][28];
		
		for(jj=0;jj<12;jj++)
		{
			n[jj] = Hills[i][jj];
		}
		
		printf("============i=%d===k[0]=%f===TNF=%f============\n",i,k[0],TNF);
		
		double delta_RIPK3 = 0.0;
		
		for(j=0;j<50;j++)
		{
			g = (1+j)*0.02;
			
			TRADD=0.0; RIPK1=0.0; Casp8=0.0; RIPK3=0.0; aCasp8=0.0;
			
			for(it=1; it<=NT; it++)
			{
				temp_TRADD=TRADD; temp_RIPK1=RIPK1; temp_Casp8=Casp8; temp_RIPK3=RIPK3; temp_aCasp8=aCasp8;
				
				t=it*dt;
				
				TRADD = temp_TRADD + dt*(k[0]*(1-temp_TRADD)*pow(TNF,n[0])/(pow(TNF,n[0])+pow(J[0],n[0])) - k[1]*temp_TRADD*pow(temp_RIPK1,n[1])/(pow(temp_RIPK1,n[1])+pow(J[1],n[1])) - d[0]*temp_TRADD);
				
				RIPK1 = temp_RIPK1 + dt*(k[2]*(g-temp_RIPK1)*pow(TNF,n[2])/(pow(TNF,n[2])+pow(J[2],n[2])) + k[3]*(g-temp_RIPK1)*pow(temp_RIPK3,n[3])/(pow(temp_RIPK3,n[3])+pow(J[3],n[3])) - k[4]*temp_RIPK1*pow(temp_TRADD,n[4])/(pow(temp_TRADD,n[4])+pow(J[4],n[4])) - k[5]*temp_RIPK1*pow(temp_Casp8,n[5])/(pow(temp_Casp8,n[5])+pow(J[5],n[5])) - d[1]*temp_RIPK1);
				
				Casp8 = temp_Casp8 + dt*(k[6]*(1-temp_Casp8-temp_aCasp8)*pow(temp_RIPK1,n[6])/(pow(temp_RIPK1,n[6])+pow(J[6],n[6])) - k[8]*temp_Casp8*pow(temp_RIPK3,n[8])/(pow(temp_RIPK3,n[8])+pow(J[8],n[8])) - d[2]*temp_Casp8);
				
				RIPK3 = temp_RIPK3 + dt*(k[9]*(1-temp_RIPK3)*pow(temp_RIPK1,n[9])/(pow(temp_RIPK1,n[9])+pow(J[9],n[9])) + k[10]*(1-temp_RIPK3)*pow(temp_RIPK3,n[10])/(pow(temp_RIPK3,n[10])+pow(J[10],n[10])) - k[11]*temp_RIPK3*pow(temp_Casp8,n[11])/(pow(temp_Casp8,n[11])+pow(J[11],n[11])) - d[3]*temp_RIPK3);
				
				aCasp8 = temp_aCasp8 + dt*(k[7]*(1-temp_Casp8-temp_aCasp8)*pow(temp_TRADD,n[7])/(pow(temp_TRADD,n[7])+pow(J[7],n[7])) - d[2]*temp_aCasp8);
			}
			
			printf("==g=%.8f=====RIPK3=%.8f\n",g,RIPK3);
			RIPK3_2h[j] = RIPK3;
			
		}
//=================choose max value================== 
		for(j=0;j<50;j++)
		{
			if(j==0)
			{
				max_RIPK3 = RIPK3_2h[0];
				maxindex = 0;
			}
			else
			{
				if(RIPK3_2h[j] > max_RIPK3)
				{
					max_RIPK3 = RIPK3_2h[j];
					maxindex = j;
				}
				
				if(j==49)
				{
					last_RIPK3 = RIPK3_2h[j];
				}
			}
		}
		//printf("%lf====%lf=====%d\n",max_RIPK3,last_RIPK3,maxindex);
//=================================================
		int sta_num=0;
		int swi = 0;
		double max_rip3, last_rip3;
		
		if(maxindex<49)
		{
			//printf("1Yes\n");
			if(maxindex==0)
			{
				//printf("2Yes\n");
				for(j=0;j<49;j++)
				{
					if(RIPK3_2h[j+1] - RIPK3_2h[j] <= 0.001)
						sta_num += 1;
				}
				//printf("%d\n",sta_num);
				if(sta_num==49 && max_RIPK3-last_RIPK3>=0.001)
				{
					printf("3Yes\n");
					
					for(j=0;j<45;j++)
					{
						if((RIPK3_2h[j+5] - RIPK3_2h[j])/0.1 >= 5)
							swi += 1;
					}
					
					if(swi > 0)
					{
						char rise_char[100];
						for(jj=0;jj<12;jj++)
						{
							sprintf(rise_char, "%.8f\t",k[jj]);
							rise_file << rise_char;
						}
						for(jj=0;jj<12;jj++)
						{
							sprintf(rise_char, "%.8f\t",J[jj]);
							rise_file << rise_char;
						}
						for(jj=0;jj<4;jj++)
						{
							sprintf(rise_char, "%.8f\t",d[jj]);
							rise_file << rise_char;
						}
						sprintf(rise_char, "%.8f\t",TNF);
						rise_file << rise_char;
						for(jj=0;jj<12;jj++)
						{
							if(jj<11)
							{
								sprintf(rise_char, "%.8f\t",n[jj]);
								rise_file << rise_char;
							}
							else
							{
								sprintf(rise_char, "%.8f",n[jj]);
								rise_file << rise_char;
							}
						}
						rise_file << std::endl;
						
						max_rip3 = 0.0;
						last_rip3 = RIPK3_2h[49];
						for(j=0;j<49;j++)
						{
							if(max_rip3 < RIPK3_2h[j])
								max_rip3 = RIPK3_2h[j];
						}
						char Hvalue_char[100];
						sprintf(Hvalue_char, "%.6f\t%.6f\t%.6f\n",max_rip3,last_rip3,max_rip3-last_rip3);
						Hvalue_file << Hvalue_char;
					}
				}
			}
			else
			{
				//printf("4Yes\n");
				for(j=0;j<maxindex;j++)
				{
					if(RIPK3_2h[j+1] - RIPK3_2h[j] >= -0.001)
						sta_num += 1;
				}
				for(j=maxindex;j<49;j++)
				{
					if(RIPK3_2h[j+1] - RIPK3_2h[j] <= 0.001)
						sta_num += 1;
				}
				//printf("%d\n",sta_num);
				if(sta_num==49 && max_RIPK3-last_RIPK3>=0.001)
				{
					printf("5Yes\n");
					
					for(j=0;j<45;j++)
					{
						if((RIPK3_2h[j+5] - RIPK3_2h[j])/0.1 >= 5)
							swi += 1;
					}
					
					if(swi > 0)
					{ 
						char rise_char[100];
						for(jj=0;jj<12;jj++)
						{
							sprintf(rise_char, "%.8f\t",k[jj]);
							rise_file << rise_char;
						}
						for(jj=0;jj<12;jj++)
						{
							sprintf(rise_char, "%.8f\t",J[jj]);
							rise_file << rise_char;
						}
						for(jj=0;jj<4;jj++)
						{
							sprintf(rise_char, "%.8f\t",d[jj]);
							rise_file << rise_char;
						}
						sprintf(rise_char, "%.8f\t",TNF);
						rise_file << rise_char;
						for(jj=0;jj<12;jj++)
						{
							if(jj<11)
							{
								sprintf(rise_char, "%.8f\t",n[jj]);
								rise_file << rise_char;
							}
							else
							{
								sprintf(rise_char, "%.8f",n[jj]);
								rise_file << rise_char;
							}
						}
						rise_file << std::endl;
						
						max_rip3 = 0.0;
						last_rip3 = RIPK3_2h[49];
						for(j=0;j<49;j++)
						{
							if(max_rip3 < RIPK3_2h[j])
								max_rip3 = RIPK3_2h[j];
						}
						char Hvalue_char[100];
						sprintf(Hvalue_char, "%.6f\t%.6f\t%.6f\n",max_rip3,last_rip3,max_rip3-last_rip3);
						Hvalue_file << Hvalue_char;
					} 
				}
			}
			
		}
		
	}//loop of NP
	
	rise_file.close();
	Hvalue_file.close();
//============================================================================================================= 
	MPI_Finalize();
	return 0;
}
