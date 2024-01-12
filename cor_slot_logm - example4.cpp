#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include<iomanip>
#include <fstream>
#include "mkl.h"
#include"SolveLinearSystem.h"

using namespace std;

#ifndef PI
#define PI 3.141592653589793
#endif

int iplus(int i, int N) {
	int ti1;
	if (i == N) {
		ti1 = 0;

	}
	else {
		ti1 = i + 1;

	}
	return ti1;
}

int iminus(int i, int N) {
	int ti2;
	if (i == 0) {

		ti2 = N;
	}
	else {
		ti2 = i - 1;
	}
	return ti2;
}

struct gridpoint
{
	double E_x;
	double E_y;

};




double slot(double e, double q, double epsl1, double epsl2,
	double rv, double vc1,double vc2, double d) {
	double dg;
	dg = -d * e + q * .5 * (1 / epsl1 - 1 / epsl2) - rv * (log(vc1) - log(vc2));
		
	return dg;
}



double Field_Plaquette_Move(struct gridpoint* grid, struct gridpoint* epsl, int N_x, int N_y,
	double delta_x, double delta_y,
	int x, int y, int t)
{
	double functional_change = .0;
	double optimal_flux_update = .0;

	double E1, E2, E3, E4 = .0;  // fields

	double w1, w2, d1, d2 = .0; // geometrical factors

	// periodic boundaries
	int x_neighbor, y_neighbor = 0;

	if (x == N_x - 1)
		x_neighbor = 0;
	else
		x_neighbor = x + 1;

	if (y == N_y - 1)
		y_neighbor = 0;
	else
		y_neighbor = y + 1;



	// Fields
	E1 = grid[x * N_y + y + t * N_x * N_y].E_x;
	E2 = grid[x_neighbor * N_y + y + t * N_x * N_y].E_y;
	E3 = grid[x * N_y + y_neighbor + t * N_x * N_y].E_x;
	E4 = grid[x * N_y + y + t * N_x * N_y].E_y;

	// geometrical factors
	d1 = delta_x;
	d2 = delta_y;



	optimal_flux_update = -delta_x * delta_y * (d1 * (E1 - E3) + d2 * (E2 - E4))
		/ (d1 * d1 * (1 / epsl[x * N_y + y].E_x + 1 / epsl[x * N_y + y_neighbor].E_x) + d2 * d2 * (1 / epsl[x * N_y + y].E_y + 1 / epsl[x_neighbor * N_y + y].E_x));
	
	if (fabs(optimal_flux_update) >= 1e-12)
	{

		functional_change = optimal_flux_update * (
			(d1 * d1 * (1 / epsl[x * N_y + y].E_x + 1 / epsl[x * N_y + y_neighbor].E_x) + d2 * d2 * (1 / epsl[x * N_y + y].E_y + 1 / epsl[x_neighbor * N_y + y].E_x)) * optimal_flux_update / (2 * delta_x * delta_y)
			+ d1 * (E1 - E3) + d2 * (E2 - E4));

		
		if (functional_change > .0) {
			printf("bad move %f\n", functional_change);
		}



		grid[x * N_y + y + t * N_x * N_y].E_x += optimal_flux_update / (d2 * epsl[x * N_y + y].E_x);
		grid[x_neighbor * N_y + y + t * N_x * N_y].E_y += optimal_flux_update / (d1 * epsl[x_neighbor * N_y + y].E_y);
		grid[x * N_y + y_neighbor + t * N_x * N_y].E_x -= optimal_flux_update / (d2 * epsl[x * N_y + y_neighbor].E_x);
		grid[x * N_y + y + t * N_x * N_y].E_y -= optimal_flux_update / (d1 * epsl[x * N_y + y].E_y);

	}


	return(functional_change);

}


double find(double a[], int n, int m)
{
	double min = a[n - 1];
	for (int i = n - 1; i < m; i++)
	{
		if (a[i] < min)
		{
			min = a[i];
			continue;
		}
	}
	return min;

}


 //log mean
double bern(double z)
{
	double B;
	if (fabs(z) > 1e-12)
		B = z / (exp(z) - 1);
	else B = 1;
	return B;
}
/*
//harmonic mean
double bern(double z)
{
	double B;
	B=2/(1+exp(z));
	return B;
}*/


double sum(double a[], int n, int m)
{
	double sum = 0;
	for (int i = n - 1; i < m; i++)
	{
		sum += a[i];
	}
	return sum;

}

double max(double a, double b)
{
	double c;
	c = fabs(a);
	if (fabs(a) < fabs(b))
	{
		c = fabs(b);
	}

	return c;

}

int main()
{
	clock_t start = clock();
	int i, j, t, k = 0;

	double l_x = 2;
	double l_y = 2;
	double T = 0.2;
	int N_x = 200;
	int N_y = N_x;
	double kappa = 0.02;
	double cbulk = 0.1;

	double delta_x = l_x / N_x;
	double delta_y = l_y / N_y;
	double delta_t =0.1*delta_x;
	double TN = T / delta_t + 1;
	int Tn = (int)(TN + 0.5);
	cout << Tn << endl;
	double q = 1 / (4.0 * 4.0 * PI * kappa * kappa);
	//double q = 0; //no born energy

	double vci1, vci2, vcj1, vcj2,vci;
	double rad1,rad2,rad0,rad;
	rad2 = 0.338; rad0 = 0.1375; rad1 = 0.358;
	double v2 = 8*rad2 * rad2 * rad2, v1 = 8 *rad1* rad1 * rad1, v0 = 8 *rad0* rad0 * rad0;
	double functional_tolrence = 1e-5;

	struct gridpoint* field = new struct gridpoint[N_x * N_y];
	double* ion1_conz = new double[N_x * N_y];
	double* ion2_conz = new double[N_x * N_y];

	struct gridpoint* J = new struct gridpoint[N_x * N_y * 2]; //Theta2
	struct gridpoint* Q = new struct gridpoint[N_x * N_y * 2];//Theta2

	double* b1 = new double[N_x * N_y];
	double* b2 = new double[N_x * N_y];
	double* rhof = new double[N_x * N_y];
	double* c1min = new double[Tn];
	double* c2min = new double[Tn];
	double* psum = new double[Tn];
	double* nsum = new double[Tn];
	double* functional = new double[Tn];
	double* PeNum1 = new double[N_x * N_y];
	double* PeNum2 = new double[N_x * N_y];
	double* Rs = new double[Tn-1];
	struct gridpoint* epsl = new struct gridpoint[N_x * N_y];

	int ti1, ti2, tj1, tj2,ti3,tj3;

	
	for (i = 0; i < N_x; i++)
	{
		ti1 = iplus(i, N_x - 1);
		for (j = 0; j < N_x; j++)
		{
			tj1 = iplus(j, N_y - 1);
			ion1_conz[i * N_y + j] = 0.1; //initial data
			ion2_conz[i * N_y + j] = 0.1;
			epsl[i * N_y + j].E_x = 1;
			epsl[i * N_y + j].E_y = 1;
			//epsl[i * N_y + j].E_x = 77 * (tanh(50 * sqrt((i * delta_x-1) * (i * delta_x-1 ) + (j * delta_y -1) * (j * delta_y -1)) - 25)
				//+ tanh(50 * sqrt((ti1 * delta_x-1) * (ti1 * delta_x-1) + (j * delta_y-1 ) * (j * delta_y-1)) - 25)) / 4 + 79 / 2;
			//epsl[i * N_y + j].E_y = 77 * (tanh(50 * sqrt((i * delta_x -1) * (i * delta_x -1) + (j * delta_y -1) * (j * delta_y-1)) - 25)
				//+ tanh(50 * sqrt((i * delta_x-1 ) * (i * delta_x-1 ) + (tj1 * delta_y-1) * (tj1 * delta_y-1)) - 25)) / 4 + 79 / 2;
			
			
		}

	}
	c1min[0] = 0.1;
	c2min[0] = 0.1;
	psum[0] = delta_x * delta_y * sum(ion1_conz, 1, N_x * N_y);
	nsum[0] = delta_x * delta_y * sum(ion2_conz, 1, N_x * N_y);

	
	double* A1 = new double[N_x * N_y * 5];
	double* A2 = new double[N_x * N_y * 5];
	double* B1 = new double[N_x * N_y * 5];
	double* B2 = new double[N_x * N_y * 5];
	double di1, di2, dj1,dj2;
	double c0;
	

	double* c1_temp = new double[N_x * N_y];
	double* c2_temp = new double[N_x * N_y];
	struct gridpoint* field_temp = new struct gridpoint[N_x * N_y*3];
	

	double theta_x, theta_y, z;
	double functional_change;

	int a1, a2, a3, a4, a5;
	a1 = 0; a2 = 1; a3 = 2; a4 = 3; a5 = 4;

	
	int r;
	MKL_INT* ia = new MKL_INT[N_x * N_y + 1];
	ia[0] = 0;
	for (i = 1; i < N_x * N_y + 1; i++) {
		ia[i] = ia[i - 1] + 5;
	}
	MKL_INT* ja = new MKL_INT[N_x * N_y * 5];
	for (i = 0; i < N_x; i++) {
		ti1 = iplus(i, N_x - 1);
		ti2 = iminus(i, N_x - 1);
		for (j = 0; j < N_y; j++) {
			tj1 = iplus(j, N_y - 1);
			tj2 = iminus(j, N_y - 1);

			if (i == 0) {
				a4 = 3; a5 = 4;
				if (j == 0) {
					a1 = 0; a2 = 1; a3 = 2;
				}
				else if (j == N_y - 1) {
					a1 = 2; a2 = 0; a3 = 1;
				}
				else {
					a1 = 1; a2 = 2; a3 = 0;
				}
			}
			else if (i == N_x - 1) {
				a4 = 0; a5 = 1;
				if (j == 0) {
					a1 = 2; a2 = 3; a3 = 4;
				}
				else if (j == N_y - 1) {
					a1 = 4; a2 = 2; a3 = 3;
				}
				else {
					a1 = 3; a2 = 4; a3 = 2;
				}
			}
			else {
				a4 = 4; a5 = 0;
				if (j == 0) {
					a1 = 1; a2 = 2; a3 = 3;
				}
				else if (j == N_y - 1) {
					a1 = 3; a2 = 1; a3 = 2;
				}
				else {
					a1 = 2; a2 = 3; a3 = 1;
				}
			}

			ja[5 * (i * N_y + j) + a1] = i * N_y + j;
			ja[5 * (i * N_y + j) + a2] = i * N_y + tj1;
			ja[5 * (i * N_y + j) + a3] = i * N_y + tj2;
			ja[5 * (i * N_y + j) + a4] = ti1 * N_y + j;
			ja[5 * (i * N_y + j) + a5] = ti2 * N_y + j;

			A1[5 * (i * N_y + j) + a1] = 2 * kappa * kappa * ((epsl[i * N_y + j].E_x + epsl[ti2 * N_y + j].E_x) / (delta_x * delta_x) + (epsl[i * N_y + j].E_y + epsl[i * N_y + tj2].E_y) / (delta_y * delta_y));
			A1[5 * (i * N_y + j) + a2] = -2 * kappa * kappa * epsl[i * N_y + j].E_y / (delta_y * delta_y);
			A1[5 * (i * N_y + j) + a3] = -2 * kappa * kappa * epsl[i * N_y + tj2].E_y / (delta_y * delta_y);
			A1[5 * (i * N_y + j) + a4] = -2 * kappa * kappa * epsl[i * N_y + j].E_x / (delta_x * delta_x);
			A1[5 * (i * N_y + j) + a5] = -2 * kappa * kappa * epsl[ti2 * N_y + j].E_x / (delta_x * delta_x);
			
			rhof[i * N_y + j] = 0;
			if (fabs((i * delta_x - 1) * (i * delta_x - 1) + (j * delta_y - 1) * (j * delta_y - 1) - 0.25) < 1e-2) {
				if (j * delta_y >= 1) {
					rhof[i * N_y + j] = 1.0;
					if (j * delta_y <= 1 && i * delta_x >= 1.5 - 1e-2)
						rhof[i * N_y + j] = -1.0;
				}

				else
					rhof[i * N_y + j] = -1.0;
			}
		}
	}
	A1[0] = 1; A1[1] = 0; A1[2] = 0; A1[3] = 0; A1[4] = 0;
	rhof[0] = 0;
	//c1_temp=initial potential
	r = sp_ds(c1_temp, ia, ja, A1, rhof, N_x * N_y);
	

	for (i = 0; i < N_x; i++)
	{
		ti1 = iplus(i, N_x - 1);
		for (j = 0; j < N_y; j++)
		{
			tj1 = iplus(j, N_y - 1);
			field[i * N_y + j].E_x = -(c1_temp[ti1 * N_y + j] - c1_temp[i * N_y + j]) / delta_x;
			field[i * N_y + j].E_y = -(c1_temp[i * N_y + tj1] - c1_temp[i * N_y + j]) / delta_y;

			field_temp[i * N_y + j].E_x = field[i * N_y + j].E_x;
			field_temp[i * N_y + j].E_y = field[i * N_y + j].E_y;
		}
	}

	functional[0] = 0;
	for (i = 0; i < N_x; i++) {
		ti1 = iplus(i, N_x - 1);
		for (j = 0; j < N_y; j++) {
			tj1 = iplus(j, N_y - 1);
			c1_temp[i * N_y + j] = ion1_conz[i * N_y + j];
			c2_temp[i * N_y + j] = ion2_conz[i * N_y + j];
			c0 = (1 - v1 * ion1_conz[i * N_y + j] - v2 * ion2_conz[i * N_y + j]) / v0;
			functional[0] += kappa * kappa * delta_x * delta_y *
				(epsl[i * N_y + j].E_x * field[i * N_y + j].E_x * field[i * N_y + j].E_x
					+ epsl[i * N_y + j].E_y * field[i * N_y + j].E_y * field[i * N_y + j].E_y)
				+ delta_x * delta_y * (ion1_conz[i * N_y + j] * log(ion1_conz[i * N_y + j])
					+ ion2_conz[i * N_y + j] * log(ion2_conz[i * N_y + j])
					+ c0 * log(c0)
					)
				+ q * delta_x * delta_y * .5 * (((ion1_conz[i * N_y + j] + ion1_conz[ti1 * N_y + j]) / rad1 + (ion2_conz[i * N_y + j] + ion2_conz[ti1 * N_y + j]) / rad2) / epsl[i * N_y + j].E_x +
					((ion1_conz[i * N_y + j] + ion1_conz[i * N_y + tj1]) / rad1 + (ion2_conz[i * N_y + j] + ion2_conz[i * N_y + tj1]) / rad2) / epsl[i * N_y + j].E_y);
			
		}
	}
	

ofstream fout("E:\\worksp\\pnp\\example1\\ex4.txt");
if (fout)
{

	for (t = 0; t < Tn - 1; t++)
	{
		functional[t + 1] = 0;
		for (i = 0; i < N_x; i++)
		{
			ti1 = iplus(i, N_x - 1); //ti1=i+1
			ti2 = iminus(i, N_x - 1);//ti2=i-1;
			ti3 = iminus(ti2, N_x - 1);//ti3=i-2;
			for (j = 0; j < N_y; j++)
			{
				tj1 = iplus(j, N_y - 1);
				tj2 = iminus(j, N_y - 1);
				tj3 = iminus(tj2, N_y - 1);
				
				if (i == 0) {
					a4 = 3; a5 = 4;
					if (j == 0) {
						a1 = 0; a2 = 1; a3 = 2;
						}
					else if (j == N_y - 1) {
						a1 = 2; a2 = 0; a3 = 1;
					}
					else {
						a1 = 1; a2 = 2; a3 = 0;
						}
				}
				else if (i == N_x - 1) {
					a4 = 0; a5 = 1;
					if (j == 0) {
						a1 = 2; a2 = 3; a3 = 4;
					}
					else if (j == N_y - 1) {
						a1 = 4; a2 = 2; a3 = 3;
					}
					else {
						a1 = 3; a2 = 4; a3 = 2;
					}
				}
				else {
					a4 = 4; a5 = 0;
					if (j == 0) {
						a1 = 1; a2 = 2; a3 = 3;
					}
					else if (j == N_y - 1) {
						a1 = 3; a2 = 1; a3 = 2;
					}
					else {
						a1 = 2; a2 = 3; a3 = 1;
					}
				}
				
				//ion1_conz
			
				vci = 1 - v1 * c1_temp[i * N_y + j] - v2 * c2_temp[i * N_y + j];
				vci1 = 1 - v1 * c1_temp[ti1 * N_y + j] - v2 * c2_temp[ti1 * N_y + j];
				vci2 = 1 - v1 * c1_temp[ti2 * N_y + j] - v2 * c2_temp[ti2 * N_y + j];
				vcj1 = 1 - v1 * c1_temp[i * N_y + tj1] - v2 * c2_temp[i * N_y + tj1];
				vcj2 = 1 - v1 * c1_temp[i * N_y + tj2] - v2 * c2_temp[i * N_y + tj2];
				z = 1.0;
				rad = rad1;
				di1 = slot(z * field_temp[i * N_y + j].E_x, q / rad, epsl[ti1 * N_y + j].E_x, epsl[ti2 * N_y + j].E_x, v1 / v0, vci1, vci, delta_x);
				di2 = slot(z * field_temp[ti2 * N_y + j].E_x, q / rad, epsl[i * N_y + j].E_x, epsl[ti3 * N_y + j].E_x, v1 / v0, vci, vci2, delta_x);
				dj1 = slot(z * field_temp[i * N_y + j].E_y, q / rad, epsl[i * N_y + tj1].E_y, epsl[i * N_y + tj2].E_y, v1 / v0, vcj1, vci, delta_y);
				dj2 = slot(z * field_temp[i * N_y + tj2].E_y, q / rad, epsl[i * N_y + j].E_y, epsl[i * N_y + tj3].E_y, v1 / v0, vci, vcj2, delta_y);
				

				A1[5 * (i * N_y + j) + a1] = 1 + kappa * delta_t * (bern(di1) + bern(-di2)) / (delta_x * delta_x)
					+ kappa * delta_t * (bern(dj1) + bern(-dj2)) / (delta_y * delta_y);
				A1[5 * (i * N_y + j) + a2] = -kappa * delta_t * (bern(-dj1)) / (delta_y * delta_y);
				A1[5 * (i * N_y + j) + a3] = -kappa * delta_t * (bern(dj2)) / (delta_y * delta_y);
				A1[5 * (i * N_y + j) + a4] = -kappa * delta_t * (bern(-di1)) / (delta_x * delta_x);
				A1[5 * (i * N_y + j) + a5] = -kappa * delta_t * (bern(di2)) / (delta_x * delta_x);
				//

				//ion2_conz
				
				z = -1.0;
				rad = rad2;

				di1 = slot(z * field_temp[i * N_y + j].E_x, q / rad, epsl[ti1 * N_y + j].E_x, epsl[ti2 * N_y + j].E_x, v2 / v0, vci1, vci, delta_x);
				di2 = slot(z * field_temp[ti2 * N_y + j].E_x, q / rad, epsl[i * N_y + j].E_x, epsl[ti3 * N_y + j].E_x, v2 / v0, vci, vci2, delta_x);
				dj1 = slot(z * field_temp[i * N_y + j].E_y, q / rad, epsl[i * N_y + tj1].E_y, epsl[i * N_y + tj2].E_y, v2 / v0, vcj1, vci, delta_y);
				dj2 = slot(z * field_temp[i * N_y + tj2].E_y, q / rad, epsl[i * N_y + j].E_y, epsl[i * N_y + tj3].E_y, v2 / v0, vci, vcj2, delta_y);

				A2[5 * (i * N_y + j) + a1] = 1 + kappa * delta_t * (bern(di1) + bern(-di2)) / (delta_x * delta_x)
					+ kappa * delta_t * (bern(dj1) + bern(-dj2)) / (delta_y * delta_y);
				A2[5 * (i * N_y + j) + a2] = -kappa * delta_t * (bern(-dj1)) / (delta_y * delta_y);
				A2[5 * (i * N_y + j) + a3] = -kappa * delta_t * (bern(dj2)) / (delta_y * delta_y);
				A2[5 * (i * N_y + j) + a4] = -kappa * delta_t * (bern(-di1)) / (delta_x * delta_x);
				A2[5 * (i * N_y + j) + a5] = -kappa * delta_t * (bern(di2)) / (delta_x * delta_x);
			}
		}
		
		r = sp_ds(ion1_conz, ia, ja, A1, c1_temp, N_x * N_y);
		
		r = sp_ds(ion2_conz, ia, ja, A2, c2_temp, N_x * N_y);

		for (i = 0; i < N_x; i++)
		{
			ti1 = iplus(i, N_x - 1);
			ti2 = iminus(i, N_x - 1);
			for (j = 0; j < N_y; j++)
			{
				tj1 = iplus(j, N_y - 1);
				tj2 = iminus(j, N_y - 1);
				c0 = (1 - v1 * ion1_conz[i * N_y + j] - v2 * ion2_conz[i * N_y + j]) / v0;
				if (ion1_conz[i * N_y + j] < 0 || ion2_conz[i * N_y + j] <0 ) {
					cout << "<0";
				}
				
				functional[t + 1] += delta_x * delta_y * ion1_conz[i * N_y + j] * log(ion1_conz[i * N_y + j])
					+ delta_x * delta_y * ion2_conz[i * N_y + j] * log(ion2_conz[i * N_y + j])
					+ delta_x * delta_y * c0 * log(c0)
					+ q * delta_x * delta_y * .5 * (((ion1_conz[i * N_y + j] + ion1_conz[ti1 * N_y + j]) / rad1 + (ion2_conz[i * N_y + j] + ion2_conz[ti1 * N_y + j]) / rad2) / epsl[i * N_y + j].E_x +
						((ion1_conz[i * N_y + j] + ion1_conz[i * N_y + tj1]) / rad1 + (ion2_conz[i * N_y + j] + ion2_conz[i * N_y + tj1]) / rad2) / epsl[i * N_y + j].E_y);
				vci = 1 - v1 * c1_temp[i * N_y + j] - v2 * c2_temp[i * N_y + j];
				vci1 = 1 - v1 * c1_temp[ti1 * N_y + j] - v2 * c2_temp[ti1 * N_y + j];
				vcj1 = 1 - v1 * c1_temp[i * N_y + tj1] - v2 * c2_temp[i * N_y + tj1];

				//ion1_conz
				
				z = 1.0; rad = rad1;

				di1 = slot(z * field_temp[i * N_y + j].E_x, q / rad, epsl[ti1 * N_y + j].E_x, epsl[ti2 * N_y + j].E_x, v1 / v0, vci1, vci, delta_x);
				dj1 = slot(z * field_temp[i * N_y + j].E_y, q / rad, epsl[i * N_y + tj1].E_y, epsl[i * N_y + tj2].E_y, v1 / v0, vcj1, vci, delta_y);
				
				
				J[i * N_y + j].E_x = -(bern(-di1)) * ion1_conz[ti1 * N_y + j] / delta_x + (bern(di1)) * ion1_conz[i * N_y + j] / delta_x;
				J[i * N_y + j].E_y = -(bern(-dj1)) * ion1_conz[i * N_y + tj1] / delta_y + (bern(dj1)) * ion1_conz[i * N_y + j] / delta_y;
				
				
				//ion2_conz
				z = -1.0; rad = rad2;
				
				
				di1 = slot(z * field_temp[i * N_y + j].E_x, q / rad, epsl[ti1 * N_y + j].E_x, epsl[ti2 * N_y + j].E_x, v2 / v0, vci1, vci, delta_x);
				dj1 = slot(z * field_temp[i * N_y + j].E_y, q / rad, epsl[i * N_y + tj1].E_y, epsl[i * N_y + tj2].E_y, v2 / v0, vcj1, vci, delta_y);

				
				Q[i * N_y + j].E_x = -(bern(-di1)) * ion2_conz[ti1 * N_y + j] / delta_x + (bern(di1)) * ion2_conz[i * N_y + j] / delta_x;
				Q[i * N_y + j].E_y = -(bern(-dj1)) * ion2_conz[i * N_y + tj1] / delta_y + (bern(dj1)) * ion2_conz[i * N_y + j] / delta_y;
				theta_x = 0; theta_y = 0;
				
				if (t > 0)
				{
					theta_x = delta_t * (J[i * N_y + j + N_x * N_y].E_x - Q[i * N_y + j + N_x * N_y].E_x) / (kappa *2* epsl[i * N_y + j].E_x) + field_temp[i * N_y + j].E_x - field_temp[i * N_y + j + N_x * N_y].E_x;
					theta_y = delta_t * (J[i * N_y + j + N_x * N_y].E_y - Q[i * N_y + j + N_x * N_y].E_y) / (kappa *2* epsl[i * N_y + j].E_y) + field_temp[i * N_y + j].E_y - field_temp[i * N_y + j + N_x * N_y].E_y;
					

				}
				field[i * N_y + j].E_x = field_temp[i * N_y + j].E_x
					+ delta_t * (-J[i * N_y + j].E_x + Q[i * N_y + j].E_x) / (kappa * 2 * epsl[i * N_y + j].E_x) + theta_x;
				field[i * N_y + j].E_y = field_temp[i * N_y + j].E_y
					+ delta_t * (-J[i * N_y + j].E_y + Q[i * N_y + j].E_y) / (kappa * 2 * epsl[i * N_y + j].E_y) + theta_y;

				functional[t + 1] += kappa * kappa * delta_x * delta_y * 
					(epsl[i * N_y + j].E_x * field[i * N_y + j].E_x * field[i * N_y + j].E_x 
						+ epsl[i * N_y + j].E_y * field[i * N_y + j].E_y * field[i * N_y + j].E_y);
                

				Q[i * N_y + j + N_x * N_y].E_x = Q[i * N_y + j].E_x;
				Q[i * N_y + j + N_x * N_y].E_y = Q[i * N_y + j].E_y;
				J[i * N_y + j + N_x * N_y].E_x = J[i * N_y + j].E_x;
				J[i * N_y + j + N_x * N_y].E_y = J[i * N_y + j].E_y;

				field_temp[i * N_y + j + N_x * N_y].E_x = field_temp[i * N_y + j].E_x;
				field_temp[i * N_y + j + N_x * N_y].E_y = field_temp[i * N_y + j].E_y;
			}

		}



		functional_change = 1.0;
		k = 0;
		while (fabs(functional_change) > functional_tolrence)
		{
			k = k + 1;
			functional_change = 0;
			for (i = 0; i < N_x; i++) {
				for (j = 0; j < N_y; j++) {
					functional_change += Field_Plaquette_Move(field, epsl, N_x, N_y, delta_x, delta_y,
						i, j, 0);
				}
			}
			functional[t + 1] += 2*kappa*kappa*functional_change;
		}
		

		for (i = 0; i < N_x; i++) {
			ti1 = iplus(i, N_x - 1);
			ti2 = iminus(i, N_x - 1);
			for (j = 0; j < N_y; j++) {
				tj1 = iplus(j, N_y - 1);
				tj2 = iminus(j, N_y - 1);
				
				vci = 1 - v1 * c1_temp[i * N_y + j] - v2 * c2_temp[i * N_y + j];
				vci1 = 1 - v1 * c1_temp[ti1 * N_y + j] - v2 * c2_temp[ti1 * N_y + j];
				vcj1 = 1 - v1 * c1_temp[i * N_y + tj1] - v2 * c2_temp[i * N_y + tj1];
				z = 1.0; rad = rad1;
				di1 = slot(z * field_temp[i * N_y + j].E_x, q / rad, epsl[ti1 * N_y + j].E_x, epsl[ti2 * N_y + j].E_x, v1 / v0, vci1, vci, delta_x);
				dj1 = slot(z * field_temp[i * N_y + j].E_y, q / rad, epsl[i * N_y + tj1].E_y, epsl[i * N_y + tj2].E_y, v1 / v0, vcj1, vci, delta_y);

				PeNum1[i * N_y + j] = max(di1,dj1);
				
				z = -1.0; rad = rad2;
				di1 = slot(z * field_temp[i * N_y + j].E_x, q / rad, epsl[ti1 * N_y + j].E_x, epsl[ti2 * N_y + j].E_x, v1 / v0, vci1, vci, delta_x);
				dj1 = slot(z * field_temp[i * N_y + j].E_y, q / rad, epsl[i * N_y + tj1].E_y, epsl[i * N_y + tj2].E_y, v1 / v0, vcj1, vci, delta_y);

				PeNum2[i * N_y + j] = max(di1,dj1);
				

				field_temp[i * N_y + j].E_x = field[i * N_y + j].E_x;
				field_temp[i * N_y + j].E_y = field[i * N_y + j].E_y;

				c1_temp[i * N_y + j] = ion1_conz[i * N_y + j];
				c2_temp[i * N_y + j] = ion2_conz[i * N_y + j];
			}
		}
		Rs[t] = k;
		c1min[t + 1] = find(ion1_conz,1, N_x * N_y);
		c2min[t + 1] = find(ion2_conz,1, N_x * N_y);
		psum[t + 1] = delta_x * delta_y * sum(ion1_conz, 1,N_x * N_y);
		nsum[t + 1] = delta_x * delta_y * sum(ion2_conz, 1, N_x * N_y);

		if (t == 0.05 / delta_t-1) {
			cout << "t == 0.05 / delta_t-1" << t << endl;
			
			for (i = 0; i < N_x; i++) {
				for (j = 0; j < N_x; j++) {
					fout << setprecision(16) << ion1_conz[i * N_y + j] << endl;
				}
			}
			
			for (i = 0; i < N_x; i++) {
				for (j = 0; j < N_x; j++) {
					fout << setprecision(16) << ion2_conz[i * N_y + j] << endl;
				}
			}
			
			for (i = 0; i < N_x; i++) {
				for (j = 0; j < N_x; j++) {
					fout << setprecision(16) << PeNum1[i * N_y + j] << endl;
				}
			}
			
			for (i = 0; i < N_x; i++) {
				for (j = 0; j < N_x; j++) {
					fout << setprecision(16) << PeNum2[i * N_y + j] << endl;
				}
			}
			for (i = 0; i < N_x; i++) {
				for (j = 0; j < N_x; j++) {
					fout << setprecision(16) <<
						sqrt(field[i * N_y + j].E_x * field[i * N_y + j].E_x * epsl[i * N_y + j].E_x * epsl[i * N_y + j].E_x
							+ field[i * N_y + j].E_y * field[i * N_y + j].E_y * epsl[i * N_y + j].E_y * epsl[i * N_y + j].E_y) << endl;
				}
			}

			
		}
		if (t == 0.1 / delta_t - 1) {
			cout << "t == 0.1 / delta_t-1" << t << endl;
			
			for (i = 0; i < N_x; i++) {
				for (j = 0; j < N_x; j++) {
					fout << setprecision(16) << ion1_conz[i * N_y + j] << endl;
				}
			}
			for (i = 0; i < N_x; i++) {
				for (j = 0; j < N_x; j++) {
					fout << setprecision(16) << ion2_conz[i * N_y + j] << endl;
				}
			}
			
			for (i = 0; i < N_x; i++) {
				for (j = 0; j < N_x; j++) {
					fout << setprecision(16) << PeNum1[i * N_y + j] << endl;
				}
			}
			
			for (i = 0; i < N_x; i++) {
				for (j = 0; j < N_x; j++) {
					fout << setprecision(16) << PeNum2[i * N_y + j] << endl;
				}
			}
			
			for (i = 0; i < N_x; i++) {
				for (j = 0; j < N_x; j++) {
					fout << setprecision(16) <<
						sqrt(field[i * N_y + j].E_x * field[i * N_y + j].E_x * epsl[i * N_y + j].E_x * epsl[i * N_y + j].E_x
							+ field[i * N_y + j].E_y * field[i * N_y + j].E_y * epsl[i * N_y + j].E_y * epsl[i * N_y + j].E_y) << endl;
				}
			}
		}

		if (t == Tn - 2) {
			
			for (i = 0; i < N_x; i++) {
				for (j = 0; j < N_x; j++) {
					fout << setprecision(16) << ion1_conz[i * N_y + j] << endl;
				}
			}
			for (i = 0; i < N_x; i++) {
				for (j = 0; j < N_x; j++) {
					fout << setprecision(16) << ion2_conz[i * N_y + j] << endl;
				}
			}
			
			for (i = 0; i < N_x; i++) {
				for (j = 0; j < N_x; j++) {
					fout << setprecision(16) << PeNum1[i * N_y + j] << endl;
				}
			}
			
			for (i = 0; i < N_x; i++) {
				for (j = 0; j < N_x; j++) {
					fout << setprecision(16) << PeNum2[i * N_y + j] << endl;
				}
			}
			for (i = 0; i < N_x; i++) {
				for (j = 0; j < N_x; j++) {
					fout << setprecision(16) <<
						sqrt(field[i * N_y + j].E_x * field[i * N_y + j].E_x * epsl[i * N_y + j].E_x * epsl[i * N_y + j].E_x
							+ field[i * N_y + j].E_y * field[i * N_y + j].E_y * epsl[i * N_y + j].E_y * epsl[i * N_y + j].E_y) << endl;
				}
			}
		}
	}


	printf("\nrunning time: %lf\n\n", 1. * (clock() - start) / CLOCKS_PER_SEC);

	

		for (t = 0; t < Tn; t++) {
			fout << setprecision(16) << c1min[t] << endl;
			fout << setprecision(16) << c2min[t] << endl;
		}

		for (t = 0; t < Tn; t++) {
			fout << setprecision(16) << psum[t] << endl;
		}
		for (t = 0; t < Tn; t++) {
			fout << setprecision(16) << functional[t] << endl;
		}
		for (t = 0; t < Tn; t++) {
			fout << setprecision(16) << nsum[t] << endl;
		}
		for (t = 0; t < Tn - 1; t++) {
			fout << Rs[t] << endl;
		}
		for (t = 0; t < Tn; t++) {
			fout << setprecision(16) << PeNum1[t] << endl;
			fout << setprecision(16) << PeNum2[t] << endl;
		}
	fout.close();
}
	delete[] field;
	delete[] ion1_conz;
	delete[] ion2_conz;
	delete[] J;
	delete[] Q;
	delete[] b1;
	delete[] b2;
	delete[] rhof;
	delete[] c1_temp;
	delete[] c2_temp;
	delete[] A1;
	delete[] A2;
	delete[] c1min;
	delete[] c2min;
	delete[] psum;
	delete[] functional;



	system("pause");
	return 0;
}

