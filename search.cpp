//2022.8.24 by xufei 
#include <fstream>
#include <iostream>
#include <typeinfo> 
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <cstring>
#include <vector>
#include<sys/file.h>

#include"mpi.h"

#define NF 3  // 方程组个数
#define NT 7.2e5
#define dt 0.001
//************************模型的参数及变量********************//
double k1, k2, k3, k4, k5, k6, k7, k8, k9;
double J1, J2, J3, J4, J5, J6, J7, J8, J9;
double n1, n2, n3, n4, n5, n6, n7, n8, n9;
double c1, c11, c2, c22, c3, c33;
double c4, c44, c5, c55, c6, c66;
double c7, c77, c8, c88, c9, c99;
double d1, d2, d3;
double x0, x1, x2;
double g, TNF; // 所有参数 

int main(int argc, char* argv[])
{
//=============================================================================================================	
	MPI_Init(&argc, &argv);//MPI库的初始化
	int numprocs, myid;//定义进程总数及进程标识
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);//获得进程数
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);//获得当前进程标识号0,1,2,3,....,numprocs - 1
//============================================================================================================= 

	int Variable;  //控制分岔参数
	
	int i, j, it, ii, jj, kk, i_k, i_c, j_c, i_np, j_s;
	
	std::ofstream output_file("Bifurcation_paras.dat",std::ofstream::out| std::ofstream::app);
	
//*==========读耦合矩阵==================
	float Couple[4697][9]={0.0};
	std::ifstream Couplefile;
	Couplefile.open("matrix.dat");
	for(i=0;i<4697;i++)
	{
		for(j=0;j<9;j++)
		{
			Couplefile >> Couple[i][j];
		}
	}
	Couplefile.close();
//===========================================*/

//*==========读二维数组文件==================	
	float paras[10000][22]={0.0};
	std::ifstream infile;
	infile.open("LHS_of_paras.dat");
	for(i=0;i<10000;i++)
	{
		for(j=0;j<22;j++)
		{
			infile >> paras[i][j];
		}
	}
	infile.close();
//===========================================*/

//*==========读二维数组文件==================	
	float Hills[10000][9]={0.0};
	std::ifstream Hillfile;
	Hillfile.open("Hill.dat");
	for(i=0;i<10000;i++)
	{
		for(j=0;j<9;j++)
		{
			Hillfile >> Hills[i][j];
		}
	}
	Hillfile.close();
	
	FILE *rise, *rise1, *rise2, *rise3, *rise4, *rise5 ;
	char fname[20];
	char fname1[20];
	char fname2[20];
	char fname3[20];
	char fname4[20];
	char fname5[20];
//===========================================*/
	for (i_c=0;i_c<1;i_c++)
	{
		
		//std::string fname = std::to_string(i_c) + "_rise.dat";
		//std::ofstream rise_file(fname.c_str(), std::ofstream::out| std::ofstream::app);
		//itoa(i, s, 10); itoa(i_c,fname,10);
		snprintf(fname, sizeof(fname), "%d", i_c);
		strcat(fname,"_rise.dat");
		rise=fopen(fname,"a");
		
		snprintf(fname1, sizeof(fname1), "%d", i_c);
		strcat(fname1,"_rise5.dat");
		rise1=fopen(fname1,"a");
		
		snprintf(fname2, sizeof(fname2), "%d", i_c);
		strcat(fname2,"_rise4-5.dat");
		rise2=fopen(fname2,"a");
		
		snprintf(fname3, sizeof(fname3), "%d", i_c);
		strcat(fname3,"_rise3-4.dat");
		rise3=fopen(fname3,"a");
		
		snprintf(fname4, sizeof(fname4), "%d", i_c);
		strcat(fname4,"_rise2-3.dat");
		rise4=fopen(fname4,"a");
		
		snprintf(fname5, sizeof(fname5), "%d", i_c);
		strcat(fname5,"_rise1-2.dat");
		rise5=fopen(fname5,"a");
		
		printf("===i_c=%d=====\n",i_c);
		
		//========000===============
		if(Couple[i_c][0]==0)
		{
			c1 = 0.0;
			c11 = 0.0;
		}
		else if(Couple[i_c][0]==1)
		{
			c1 = 1.0;
			c11 = 0.0;
		}
		else
		{
			c1 = 0.0;
			c11 = -1.0;
		}
		//==========================
		//========111===============
		if(Couple[i_c][1]==0)
		{
			c2 = 0.0;
			c22 = 0.0;
		}
		else if(Couple[i_c][1]==1)
		{
			c2 = 1.0;
			c22 = 0.0;
		}
		else
		{
			c2 = 0.0;
			c22 = -1.0;
		}
		//==========================
		//========222===============
		if(Couple[i_c][2]==0)
		{
			c3 = 0.0;
			c33 = 0.0;
		}
		else if(Couple[i_c][2]==1)
		{
			c3 = 1.0;
			c33 = 0.0;
		}
		else
		{
			c3 = 0.0;
			c33 = -1.0;
		}
		//==========================
		//========333===============
		if(Couple[i_c][3]==0)
		{
			c4 = 0.0;
			c44 = 0.0;
		}
		else if(Couple[i_c][3]==1)
		{
			c4 = 1.0;
			c44 = 0.0;
		}
		else
		{
			c4 = 0.0;
			c44 = -1.0;
		}
		//==========================
		//========444===============
		if(Couple[i_c][4]==0)
		{
			c5 = 0.0;
			c55 = 0.0;
		}
		else if(Couple[i_c][4]==1)
		{
			c5 = 1.0;
			c55 = 0.0;
		}
		else
		{
			c5 = 0.0;
			c55 = -1.0;
		}
		//==========================
		//========555===============
		if(Couple[i_c][5]==0)
		{
			c6 = 0.0;
			c66 = 0.0;
		}
		else if(Couple[i_c][5]==1)
		{
			c6 = 1.0;
			c66 = 0.0;
		}
		else
		{
			c6 = 0.0;
			c66 = -1.0;
		}
		//==========================
		//========666===============
		if(Couple[i_c][6]==0)
		{
			c7 = 0.0;
			c77 = 0.0;
		}
		else if(Couple[i_c][6]==1)
		{
			c7 = 1.0;
			c77 = 0.0;
		}
		else
		{
			c7 = 0.0;
			c77 = -1.0;
		}
		//==========================
		//========777===============
		if(Couple[i_c][7]==0)
		{
			c8 = 0.0;
			c88 = 0.0;
		}
		else if(Couple[i_c][7]==1)
		{
			c8 = 1.0;
			c88 = 0.0;
		}
		else
		{
			c8 = 0.0;
			c88 = -1.0;
		}
		//==========================
		//========888===============
		if(Couple[i_c][8]==0)
		{
			c9 = 0.0;
			c99 = 0.0;
		}
		else if(Couple[i_c][8]==1)
		{
			c9 = 1.0;
			c99 = 0.0;
		}
		else
		{
			c9 = 0.0;
			c99 = -1.0;
		}
		//==========================
		double num_biphasic = 0.0;
		double tot_biphasic = 0.0;
		
		for (kk = myid; kk < 5000; kk += numprocs)
		{
			printf("===kk=%d=====\n",kk);
			
			k1=paras[kk][0];
			k2=paras[kk][1];
			k3=paras[kk][2];
			k4=paras[kk][3];
			k5=paras[kk][4];
			k6=paras[kk][5];
			k7=paras[kk][6];
			k8=paras[kk][7];
			k9=paras[kk][8];
			
			J1=paras[kk][9];
			J2=paras[kk][10];
			J3=paras[kk][11];
			J4=paras[kk][12];
			J5=paras[kk][13];
			J6=paras[kk][14];
			J7=paras[kk][15];
			J8=paras[kk][16];
			J9=paras[kk][17];
			
			d1=paras[kk][18];
			d2=paras[kk][19];
			d3=paras[kk][20];
			
			n1 = Hills[kk][0];
			n2 = Hills[kk][1];
			n3 = Hills[kk][2];
			n4 = Hills[kk][3];
			n5 = Hills[kk][4];
			n6 = Hills[kk][5];
			n7 = Hills[kk][6];
			n8 = Hills[kk][7];
			n9 = Hills[kk][8];
			
			TNF=paras[kk][21]; 
			
			//printf("===kk=%d==k1=%lf==J1=%lf===\n",kk,k1,J1);
			
			std::vector<double> rip(20,0.0);
			
			double max_value = 0.0;
			double last_value = 0.0;
			int maxindex = 0;
			
			for (Variable = 0; Variable < 20; Variable++)
			{
				g=0.05*(Variable+1); //分岔参数
				
				double x[NF] = {0.0, 0.0, 0.0}; // 变量初值
				
				std::vector<double> y(1800,0.0);
				
				jj = -1;
				
				for(it=0;it<NT;it++)
				{
					x0 = c1*k1*(g-x[0])*pow(x[0],n1)/(pow(x[0],n1)+pow(J1,n1)) + c11*k1*x[0]*pow(x[0],n1)/(pow(x[0],n1)+pow(J1,n1)) + c2*k2*(g-x[0])*pow(x[1],n2)/(pow(x[1],n2)+pow(J2,n2)) + c22*k2*x[0]*pow(x[1],n2)/(pow(x[1],n2)+pow(J2,n2)) + c3*k3*(g-x[0])*pow(x[2],n3)/(pow(x[2],n3)+pow(J3,n3)) + c33*k3*x[0]*pow(x[2],n3)/(pow(x[2],n3)+pow(J3,n3)) + TNF*(g-x[0]) - d1*x[0];
					x1 = c4*k4*(1-x[1])*pow(x[0],n4)/(pow(x[0],n4)+pow(J4,n4)) + c44*k4*x[1]*pow(x[0],n4)/(pow(x[0],n4)+pow(J4,n4)) + c5*k5*(1-x[1])*pow(x[1],n5)/(pow(x[1],n5)+pow(J5,n5)) + c55*k5*x[1]*pow(x[1],n5)/(pow(x[1],n5)+pow(J5,n5)) + c6*k6*(1-x[1])*pow(x[2],n6)/(pow(x[2],n6)+pow(J6,n6)) + c66*k6*x[1]*pow(x[2],n6)/(pow(x[2],n6)+pow(J6,n6)) - d2*x[1];
					x2 = c7*k7*(1-x[2])*pow(x[0],n7)/(pow(x[0],n7)+pow(J7,n7)) + c77*k7*x[2]*pow(x[0],n7)/(pow(x[0],n7)+pow(J7,n7)) + c8*k8*(1-x[2])*pow(x[1],n8)/(pow(x[1],n8)+pow(J8,n8)) + c88*k8*x[2]*pow(x[1],n8)/(pow(x[1],n8)+pow(J8,n8)) + c9*k9*(1-x[2])*pow(x[2],n9)/(pow(x[2],n9)+pow(J9,n9)) + c99*k9*x[2]*pow(x[2],n9)/(pow(x[2],n9)+pow(J9,n9)) - d3*x[2];
					
					x[0] = x[0] + x0*dt;
					x[1] = x[1] + x1*dt;
					x[2] = x[2] + x2*dt;
					
					if(it>=5.4e5 && it%100==0)
					{
						jj += 1;
						y[jj] = x[2];
					}
				}
				int num_peak=0;
				
				for(ii=1;ii<1799;ii++)
				{
					if(y[ii]-y[ii-1]>0.0001 && y[ii]-y[ii+1]>0.0001)
						num_peak += 1;
				}
				
				if(num_peak<5)
				{
					rip[Variable] = x[2];
				}
				else
				{
					for(ii=0;ii<1800;ii++)
					{
						rip[Variable] += y[ii];
					}
					rip[Variable] = rip[Variable]/1800.0;
				}
				
			}
//==============================================================================================			
			for(ii=0;ii<20;ii++)
			{
				if(ii==0)
				{
					max_value = rip[ii];
					maxindex = 0;
				}
				else
				{
					if(rip[ii] > max_value)
					{
						max_value = rip[ii];
						maxindex = ii;
					}
					
					if(ii==19)
					{
						last_value = rip[ii];
					}
				}
			}
//==============================================================================================
			int sta_num=0;
			int swi = 0;
			int swi4 = 0;
			int swi3 = 0;
			int swi2 = 0;
			int swi1 = 0;
			if(maxindex<19)
			{
				//printf("1Yes\n");
				if(maxindex==0)
				{
					
					for(ii=0;ii<19;ii++)
					{
						if(rip[ii+1] - rip[ii] <= 0.001)
							sta_num += 1;
					}
					
					if(sta_num==19 && max_value-last_value>=0.01)
					{
						printf("3Yes\n");
						if(0 == flock(fileno(rise),LOCK_EX))
						{
							for(i_k=0;i_k<22;i_k++)
							{
								fprintf(rise,"%.6lf\t",paras[kk][i_k]);
							}
							for(i_k=0;i_k<9;i_k++)
							{
								if(i_k<8)
									fprintf(rise,"%.6lf\t",Hills[kk][i_k]);
								else
									fprintf(rise,"%.6lf",Hills[kk][i_k]);
							}
							fprintf(rise,"\n");
							
							flock(fileno(rise),LOCK_UN);
						}
//======================================================================================= 
						for(j_s=0;j_s<18;j_s++)
						{
							if((rip[j_s+2] - rip[j_s])/0.1 >= 5)
								swi += 1;
						}
						
						if(swi > 0)
						{
							num_biphasic += 1.0;
							
							if(0 == flock(fileno(rise1),LOCK_EX))
							{
								for(i_k=0;i_k<22;i_k++)
								{
									fprintf(rise1,"%.6lf\t",paras[kk][i_k]);
								}
								for(i_k=0;i_k<9;i_k++)
								{
									if(i_k<8)
										fprintf(rise1,"%.6lf\t",Hills[kk][i_k]);
									else
										fprintf(rise1,"%.6lf",Hills[kk][i_k]);
								}
								fprintf(rise1,"\n");
								
								flock(fileno(rise1),LOCK_UN);
							}
						}
//======================================================================================= 

//======================================================================================= 
						for(j_s=0;j_s<18;j_s++)
						{
							if((rip[j_s+2] - rip[j_s])/0.1 >= 4 && (rip[j_s+2] - rip[j_s])/0.1 < 5)
								swi4 += 1;
						}
						
						if(swi4 > 0)
						{
							if(0 == flock(fileno(rise2),LOCK_EX))
							{
								for(i_k=0;i_k<22;i_k++)
								{
									fprintf(rise2,"%.6lf\t",paras[kk][i_k]);
								}
								for(i_k=0;i_k<9;i_k++)
								{
									if(i_k<8)
										fprintf(rise2,"%.6lf\t",Hills[kk][i_k]);
									else
										fprintf(rise2,"%.6lf",Hills[kk][i_k]);
								}
								fprintf(rise2,"\n");
								
								flock(fileno(rise2),LOCK_UN);
							}
						}
//=======================================================================================

//======================================================================================= 
						for(j_s=0;j_s<18;j_s++)
						{
							if((rip[j_s+2] - rip[j_s])/0.1 >= 3 && (rip[j_s+2] - rip[j_s])/0.1 < 4)
								swi3 += 1;
						}
						
						if(swi3 > 0)
						{
							if(0 == flock(fileno(rise3),LOCK_EX))
							{
								for(i_k=0;i_k<22;i_k++)
								{
									fprintf(rise3,"%.6lf\t",paras[kk][i_k]);
								}
								for(i_k=0;i_k<9;i_k++)
								{
									if(i_k<8)
										fprintf(rise3,"%.6lf\t",Hills[kk][i_k]);
									else
										fprintf(rise3,"%.6lf",Hills[kk][i_k]);
								}
								fprintf(rise3,"\n");
								
								flock(fileno(rise3),LOCK_UN);
							}
						}
//=======================================================================================
//======================================================================================= 
						for(j_s=0;j_s<18;j_s++)
						{
							if((rip[j_s+2] - rip[j_s])/0.1 >= 2 && (rip[j_s+2] - rip[j_s])/0.1 < 3)
								swi2 += 1;
						}
						
						if(swi2 > 0)
						{
							if(0 == flock(fileno(rise4),LOCK_EX))
							{
								for(i_k=0;i_k<22;i_k++)
								{
									fprintf(rise4,"%.6lf\t",paras[kk][i_k]);
								}
								for(i_k=0;i_k<9;i_k++)
								{
									if(i_k<8)
										fprintf(rise4,"%.6lf\t",Hills[kk][i_k]);
									else
										fprintf(rise4,"%.6lf",Hills[kk][i_k]);
								}
								fprintf(rise4,"\n");
								
								flock(fileno(rise4),LOCK_UN);
							}
						}
//=======================================================================================
//======================================================================================= 
						for(j_s=0;j_s<18;j_s++)
						{
							if((rip[j_s+2] - rip[j_s])/0.1 >= 1 && (rip[j_s+2] - rip[j_s])/0.1 < 2)
								swi1 += 1;
						}
						
						if(swi1 > 0)
						{
							if(0 == flock(fileno(rise5),LOCK_EX))
							{
								for(i_k=0;i_k<22;i_k++)
								{
									fprintf(rise5,"%.6lf\t",paras[kk][i_k]);
								}
								for(i_k=0;i_k<9;i_k++)
								{
									if(i_k<8)
										fprintf(rise5,"%.6lf\t",Hills[kk][i_k]);
									else
										fprintf(rise5,"%.6lf",Hills[kk][i_k]);
								}
								fprintf(rise5,"\n");
								
								flock(fileno(rise5),LOCK_UN);
							}
						}
//=======================================================================================

					}
				}
				else
				{
					for(ii=0;ii<maxindex;ii++)
					{
						if(rip[ii+1] - rip[ii] >= -0.001)
							sta_num += 1;
					}
					for(ii=maxindex;ii<19;ii++)
					{
						if(rip[ii+1] - rip[ii] <= 0.001)
							sta_num += 1;
					}
					
					if(sta_num==19 && max_value-last_value>=0.01)
					{
						printf("5Yes\n");
						
						if(0 == flock(fileno(rise),LOCK_EX))
						{
							for(i_k=0;i_k<22;i_k++)
							{
								fprintf(rise,"%.6lf\t",paras[kk][i_k]);
							}
							for(i_k=0;i_k<9;i_k++)
							{
								if(i_k<8)
									fprintf(rise,"%.6lf\t",Hills[kk][i_k]);
								else
									fprintf(rise,"%.6lf",Hills[kk][i_k]);
							}
							fprintf(rise,"\n");
							
							flock(fileno(rise5),LOCK_UN);
						}
//=======================================================================================	
						for(j_s=0;j_s<18;j_s++)
						{
							if((rip[j_s+2] - rip[j_s])/0.1 >= 5)
								swi += 1;
						}
						
						if(swi > 0)
						{
							num_biphasic += 1.0;
							
							if(0 == flock(fileno(rise1),LOCK_EX))
							{
								for(i_k=0;i_k<22;i_k++)
								{
									fprintf(rise1,"%.6lf\t",paras[kk][i_k]);
								}
								for(i_k=0;i_k<9;i_k++)
								{
									if(i_k<8)
										fprintf(rise1,"%.6lf\t",Hills[kk][i_k]);
									else
										fprintf(rise1,"%.6lf",Hills[kk][i_k]);
								}
								fprintf(rise1,"\n");
								
								flock(fileno(rise1),LOCK_UN);
							}
						}
//======================================================================================= 
						for(j_s=0;j_s<18;j_s++)
						{
							if((rip[j_s+2] - rip[j_s])/0.1 >= 4 && (rip[j_s+2] - rip[j_s])/0.1 < 5)
								swi4 += 1;
						}
						
						if(swi4 > 0)
						{
							if(0 == flock(fileno(rise2),LOCK_EX))
							{
								for(i_k=0;i_k<22;i_k++)
								{
									fprintf(rise2,"%.6lf\t",paras[kk][i_k]);
								}
								for(i_k=0;i_k<9;i_k++)
								{
									if(i_k<8)
										fprintf(rise2,"%.6lf\t",Hills[kk][i_k]);
									else
										fprintf(rise2,"%.6lf",Hills[kk][i_k]);
								}
								fprintf(rise2,"\n");
								
								flock(fileno(rise2),LOCK_UN);
							}
						}
//=======================================================================================
//======================================================================================= 
						for(j_s=0;j_s<18;j_s++)
						{
							if((rip[j_s+2] - rip[j_s])/0.1 >= 3 && (rip[j_s+2] - rip[j_s])/0.1 < 4)
								swi3 += 1;
						}
						
						if(swi3 > 0)
						{
							if(0 == flock(fileno(rise3),LOCK_EX))
							{
								for(i_k=0;i_k<22;i_k++)
								{
									fprintf(rise3,"%.6lf\t",paras[kk][i_k]);
								}
								for(i_k=0;i_k<9;i_k++)
								{
									if(i_k<8)
										fprintf(rise3,"%.6lf\t",Hills[kk][i_k]);
									else
										fprintf(rise3,"%.6lf",Hills[kk][i_k]);
								}
								fprintf(rise3,"\n");
								
								flock(fileno(rise3),LOCK_UN);
							}
						}
//=======================================================================================
//======================================================================================= 
						for(j_s=0;j_s<18;j_s++)
						{
							if((rip[j_s+2] - rip[j_s])/0.1 >= 2 && (rip[j_s+2] - rip[j_s])/0.1 < 3)
								swi2 += 1;
						}
						
						if(swi2 > 0)
						{
							if(0 == flock(fileno(rise4),LOCK_EX))
							{
								for(i_k=0;i_k<22;i_k++)
								{
									fprintf(rise4,"%.6lf\t",paras[kk][i_k]);
								}
								for(i_k=0;i_k<9;i_k++)
								{
									if(i_k<8)
										fprintf(rise4,"%.6lf\t",Hills[kk][i_k]);
									else
										fprintf(rise4,"%.6lf",Hills[kk][i_k]);
								}
								fprintf(rise4,"\n");
								
								flock(fileno(rise4),LOCK_UN);
							}
						}
//=======================================================================================
//======================================================================================= 
						for(j_s=0;j_s<18;j_s++)
						{
							if((rip[j_s+2] - rip[j_s])/0.1 >= 1 && (rip[j_s+2] - rip[j_s])/0.1 < 2)
								swi1 += 1;
						}
						
						if(swi1 > 0)
						{
							if(0 == flock(fileno(rise5),LOCK_EX))
							{
								for(i_k=0;i_k<22;i_k++)
								{
									fprintf(rise5,"%.6lf\t",paras[kk][i_k]);
								}
								for(i_k=0;i_k<9;i_k++)
								{
									if(i_k<8)
										fprintf(rise5,"%.6lf\t",Hills[kk][i_k]);
									else
										fprintf(rise5,"%.6lf",Hills[kk][i_k]);
								}
								fprintf(rise5,"\n");
								
								flock(fileno(rise5),LOCK_UN);
							}
						}
//=======================================================================================
					}
				}
				
			}
		}
		fclose(rise);
		fclose(rise1);
		fclose(rise2);
		fclose(rise3);
		fclose(rise4);
		fclose(rise5);
		
		MPI_Reduce(&num_biphasic, &tot_biphasic, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		if(myid == 0)
		{
			printf("=====Probability=%f==============\n", tot_biphasic/5000.0);
			
			char out_char[100];
			sprintf(out_char, "%f\t",tot_biphasic/5000.0);
			output_file << out_char;
			for (j = 0; j < 9; j++)
			{
				sprintf(out_char, "%f\t", Couple[i_c][j]);//要比较数值的小数点位数，同上
				output_file << out_char;
			}
			output_file << std::endl;
		}
	}
	output_file.close();
//============================================================================================================= 
	MPI_Finalize();
	return 0;
}
