#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "winbgi2.h"
#include "solver.h"

const unsigned int mx1 = 40, my1 = 10, mx2 = 10, my2 = 20; //domain dimensions
const unsigned int n = mx1 * (my1 + my2); //liczba wezlow
enum boundary_type { INTERIOR, LEFT, RIGHT, UP, DOWN, CORNER_LD, CORNER_RD, CORNER_LU, CORNER_RU, EMPTY }; //types of boundary conditions
enum boundary_type boundary[n]; //boundary conditions
const double D = 3.92; //diffusivity [m^2/godz]
const double beta = 30; //heat conductivity [W/(m*K)]
const double h = 1; //mesh step
const double T = 200; //number of time steps
const double dt = 0.1; //time step
const double C0 = 90; //control starting value
const double V = 80; //target temperature
const int max_iter = 100;
const double max_temp = 200;

int node(int x, int y)
{
	return y * mx1 + x;
}

void boundary_conditions(double** A)
{
	for (int i = 0; i < n; i++)
		boundary[i] = EMPTY;

	//lower beam
	for (int i = 1; i < mx1 - 1; i++)
		for (int j = 1; j < my1 - 1; j++)
			boundary[node(i, j)] = INTERIOR;
	boundary[node(0, 0)] = CORNER_LD;
	boundary[node(mx1 - 1, 0)] = CORNER_RD;
	boundary[node(0, my1 - 1)] = CORNER_LU;
	boundary[node(mx1 - 1, my1 - 1)] = CORNER_RU;
	for (int j = 1; j < my1-1; j++)
	{
		boundary[node(0, j)] = LEFT;
		boundary[node(mx1 - 1, j)] = RIGHT;
	}
	for (int i = 1; i < mx1-1; i++)
		boundary[node(i, 0)] = DOWN;
	for (int i = 1; i < (mx1 - mx2) / 2; i++)
		boundary[node(i, my1 - 1)] = UP;
	for (int i = (mx1 - mx2) / 2; i < mx1 - 1; i++)
		boundary[node(i, my1 - 1)] = UP;

	//upper beam
	for (int i = (mx1 - mx2) / 2; i < (mx1 + mx2) / 2; i++)
		for (int j = my1 - 1; j < my1 + my2 - 1; j++)
			boundary[node(i, j)] = INTERIOR;
	boundary[node((mx1 - mx2) / 2, my1 + my2 - 1)] = CORNER_LU;
	boundary[node((mx1 + mx2) / 2, my1 + my2 - 1)] = CORNER_RU;
	for (int j = my1; j < my1 + my2 - 1; j++)
	{
		boundary[node((mx1 - mx2) / 2, j)] = LEFT;
		boundary[node((mx1 + mx2) / 2, j)] = RIGHT;
	}
	for (int i = (mx1 - mx2) / 2 + 1; i < (mx1 + mx2) / 2; i++)
		boundary[node(i, my1 + my2 - 1)] = UP;
}

void fill_matrix(double** A)
{
	double alfa = dt * D / (h * h);
	for (int i = 0; i < mx1; i++)
	{
		for (int j = 0; j < my1 + my2; j++)
		{
			if (boundary[node(i, j)] == INTERIOR)
			{
				A[node(i, j)][node(i, j)] = 1 + 4 * alfa;
				A[node(i, j)][node(i - 1, j)] = -alfa;
				A[node(i, j)][node(i + 1, j)] = -alfa;
				A[node(i, j)][node(i, j - 1)] = -alfa;
				A[node(i, j)][node(i, j + 1)] = -alfa;
			}
			else if (boundary[node(i, j)] == LEFT)
			{
				A[node(i, j)][node(i, j)] = 1 + 4 * alfa;
				A[node(i, j)][node(i + 1, j)] = -2 * alfa;
				A[node(i, j)][node(i, j - 1)] = -alfa;
				A[node(i, j)][node(i, j + 1)] = -alfa;
			}
			else if (boundary[node(i, j)] == RIGHT)
			{
				A[node(i, j)][node(i, j)] = 1 + 4 * alfa;
				A[node(i, j)][node(i - 1, j)] = -2 * alfa;
				A[node(i, j)][node(i, j - 1)] = -alfa;
				A[node(i, j)][node(i, j + 1)] = -alfa;
			}
			else if (boundary[node(i, j)] == UP)
			{
				A[node(i, j)][node(i, j)] = 1 + 4 * alfa;
				A[node(i, j)][node(i - 1, j)] = -alfa;
				A[node(i, j)][node(i + 1, j)] = -alfa;
				A[node(i, j)][node(i, j - 1)] = -2 * alfa;
			}
			else if (boundary[node(i, j)] == DOWN)
			{
				A[node(i, j)][node(i, j)] = 1 + 4 * alfa;
				A[node(i, j)][node(i - 1, j)] = -alfa;
				A[node(i, j)][node(i + 1, j)] = -alfa;
				A[node(i, j)][node(i, j + 1)] = -2 * alfa;
			}
			else if (boundary[node(i, j)] == CORNER_LD)
			{
				A[node(i, j)][node(i, j)] = 1 + 4 * alfa;
				A[node(i, j)][node(i + 1, j)] = -2 * alfa;
				A[node(i, j)][node(i, j + 1)] = -2 * alfa;
			}
			else if (boundary[node(i, j)] == CORNER_RD)
			{
				A[node(i, j)][node(i, j)] = 1 + 4 * alfa;
				A[node(i, j)][node(i - 1, j)] = -2 * alfa;
				A[node(i, j)][node(i, j + 1)] = -2 * alfa;
			}
			else if (boundary[node(i, j)] == CORNER_LU)
			{
				A[node(i, j)][node(i, j)] = 1 + 4 * alfa;
				A[node(i, j)][node(i + 1, j)] = -2 * alfa;
				A[node(i, j)][node(i, j - 1)] = -2 * alfa;
			}
			else if (boundary[node(i, j)] == CORNER_RU)
			{
				A[node(i, j)][node(i, j)] = 1 + 4 * alfa;
				A[node(i, j)][node(i - 1, j)] = -2 * alfa;
				A[node(i, j)][node(i, j - 1)] = -2 * alfa;
			}
			else
				A[node(i, j)][node(i, j)] = 1;
		}
	}
}

void quad_color(int x, int y, double temp, double max, double min)
{
	double px[4], py[4];
	px[0] = x;
	py[0] = y;
	px[1] = x + 1;
	py[1] = y;
	px[2] = x + 1;
	py[2] = y + 1;
	px[3] = x;
	py[3] = y + 1;
	const double color_scaled = (temp - min) / (max - min);
	setcolor(color_scaled);
	polygon(px, py, 4);
	setgray(1);
	line(px[0], py[0], px[1], py[1]);
	line(px[1], py[1], px[2], py[2]);
	line(px[2], py[2], px[3], py[3]);
	line(px[3], py[3], px[0], py[0]);
}

void draw(double* x)
{
	double max = 90;
	double min = 15;

	for (int i = 0; i < mx1; i++)
		for (int j = 0; j < my1 + my2; j++)
			if (boundary[node(i, j)] != EMPTY)
				quad_color(i, j, x[node(i, j)], max, min);
}

void primal(double**A, double*x, double* C)
{
	double* b = (double*)malloc(n * sizeof(double));
	double alfa = dt * D / (h * h);
	for (int i = 0; i < n; i++)
		x[i] = 15;
	for (int t = 0; t < T; t++)
	{
		for (int i = 0; i < mx1; i++)
		{
			for (int j = 0; j < my1 + my2; j++)
			{
				b[node(i, j)] = x[node(i, j)];
				if(j == 0 || j == my1 + my2 - 1)
					b[node(i,j)] += dt * alfa * 2 * h * beta / D * (C[t] - x[node(i,j)]);
			}
		}
		solver_CG(n, A, b, x, x);

		draw(x);
		animate(10);
	}
	free(b);
}

void adjoint(double** A, double* u, double* C)
{
	double* b = (double*)malloc(n * sizeof(double));
	double step = 0.1;
	double alfa = dt * D / (h * h);

	for (int t = T - 1; t >= 0; t--)
	{
		for (int i = 0; i < mx1; i++)
		{
			for (int j = 0; j < my1 + my2; j++)
			{
				b[node(i, j)] = u[node(i, j)];
				if (j == 0 || j == my1 + my2 - 1)
					b[node(i, j)] -= dt * alfa * beta * u[node(i, j)];
			}
		}

		solver_CG(n, A, b, u, u);
		int count = 0;
		double norm_u = 0;
		for (int i = 0; i < n; i++)
			if (boundary[i] != EMPTY)
			{
				norm_u += u[i];
				count++;
			}
		norm_u /= count;
		C[t] += step * beta * norm_u;
		if (C[t] > max_temp)
			C[t] = max_temp;
	}
	free(b);
}

int main()
{
	//allocation of variables and matrixes
	double** A = (double**)malloc(n * sizeof(double*));
	for (int i = 0; i < n; i++)
		A[i] = (double*)malloc(n * sizeof(double));
	double* x = (double*)malloc(n * sizeof(double));
	double* u = (double*)malloc(n * sizeof(double)); 
	double* C = (double*)malloc(T * sizeof(double));

	double alfa = dt * D / (h * h);

	//adding initial conditions
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			A[i][j] = 0;
		x[i] = 15;
	}
	for (int t = 0; t < T; t++)
		C[t] = C0;

	//adding boundary conditions
	boundary_conditions(A);

	//filling matrix A
	fill_matrix(A);

	graphics(600, 600);
	scale(-5, -10, mx1 + 5, my1 + my2 + 10);

	double res_u = 2; //temperature difference residual
	int iter = 0;
	int count = 0;

	while ((iter < max_iter && res_u > 1))
	{
		primal(A, x, C); //solving primal problem

		for (int i = 0; i < n; i++)
			u[i] = V - x[i]; 

		//printing temp diff residual
		res_u = 0;
		for (int i = 0; i < n; i++)
			if (boundary[i] != EMPTY)
			{
				res_u += u[i]*u[i];
			}
		res_u = sqrt(res_u);
		printf("temperature difference residual after %d iteration is: %lf\n", iter+1, res_u);
		adjoint(A, u, C); //solving adjoint problem
		iter++;
	}

	//plotting control temperature
	wait();
	clear();
	scale(0, 0, T, max_temp);
	for (int t = 0; t < T-1; t++)
	{
		line(t, C[t], t + 1, C[t + 1]);
	}
	
	wait(); wait();
	//freeing memory
	for (int i = 0; i < n; i++)
		free(A[i]);
	free(A); free(x); free(u); free(C);
	return 0;
}