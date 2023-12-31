/* This is a sample code on how to simulate lattice 
 * boltzmann method on a simple 2D lattice. 
 *
 * In Lattice boltzmann methods we simulate the 
 * time evolution of one-particle distribution - 
 * which can be seen as the average number of particles 
 * in a region. In lattice boltzmann methods, the system 
 * is divided into a lattice, where each site represents a 
 * location in the system. Along with this coarse graining 
 * of the positions, the velocities are also discretized,
 * depending on the lattice connectivity. 
 *
 * If a site in the lattice is connected to n other 
 * sites, then we have n velocities. 
 *
 * varghese babu 
 *
 * 17/12/2023 
 *
 */
#include<iostream>
#include<math.h>

int main()
{
	// these are the sizes along x direction and y direction
	int size_x = 20; 
	int size_y = 10; 

	float m=1.;
	float k=1;	
	float T=1.;
	float c = sqrt(3.*k*T/m); // this is like the average velocity (check??!!)

	float w = 1.; // this is related to updating of the distribution in the collision step

	// lattice connectivity 
	int n_connect = 9;

	float * densities = nullptr;
	
	densities = new float[n_connect*size_x*size_y];

	float * eq_densities = nullptr;
	
	eq_densities = new float[n_connect*size_x*size_y];

	float * eq_densities = nullptr;
	temp_densities = new float[n_connect*size_x*size_y]; // I need these. 
							  
	int * solid = nullptr;

	solid = new int[size_x*size_y];
	// initialization
	
	for(int x=0; x<size_x; x++)
	{
		for(int y=0; y<size_y; y++)
		{
			if(y==size_y-1 or y==0)
			{
				solid[y*size_x+x]=1.;
			}
			else 
				solid[y*size_x+x]=0.;
		}
	}
	for(int n=0; n<n_connect; n++)
	{
		for(int x=0; x<size_x; x++)
		{
			for(int y=0; y<size_y; y++)
			{
				densities[n*size_x*size_y+y*size_x+x]=1.;
			}
		}
	}

	// n_connect is the number of velocities in the system. 
	
	// the following are the vectors which corresponds to the 
	// n_connect velocities. 

	int * e_vec_x = new int[n_connect];
	int * e_vec_y = new int[n_connect];
	
	int n_opposite = new int[n_connect];

	n_opposite[0]=0;
	n_opposite[1]=3;
	n_opposite[2]=4;
	n_opposite[3]=1;
	n_opposite[4]=2;
	n_opposite[5]=8;
	n_opposite[6]=7;
	n_opposite[7]=6;
	n_opposite[8]=5;

	e_vec_x[0] = 0;
	e_vec_y[0] = 0;  

	e_vec_x[1] = 1;
	e_vec_y[1] = 0; 

	e_vec_x[2] = 0;
	e_vec_y[2] = 1; // this is flowing up 

	e_vec_x[3] = -1;
	e_vec_y[3] = 0; 

	e_vec_x[4] = 0;
	e_vec_y[4] = -1;  // this is flowing down

	e_vec_x[5] = 1;
	e_vec_y[5] = 1; // this is flowing up

	e_vec_x[6] = 1;
	e_vec_y[6] = -1; // this is flowing down

	e_vec_x[7] = -1;
	e_vec_y[7] = 1; // this is flowing up

	e_vec_x[8] = -1;
	e_vec_y[8] = -1; // this is flowing down

	// we also have the weights associated with these vectors. 

	float * weights = new float[n_connect];

	weights[0]=4./9.;

	weights[1]=1./9.;
	weights[2]=1./9.;
	weights[3]=1./9.;
	weights[4]=1./9.;

	weights[5]=1./36.;
	weights[6]=1./36.;
	weights[7]=1./36.;
	weights[8]=1./36.;

	// we need a few more arrays to store the average velocity 
	// and the density

	float * ux=new float[size_x*size_y];
	float * uy=new float[size_x*size_y];
	float * rho=new float[size_x*size_y];
	// now we do the simulation 

	int timesteps = 10;

	for(int t=0; t<timesteps;t++)
	{
		// first is the streaming step ??.
		for(int x=0; x<size_x; x++)
		{
			for(int y=0; y<size_y; y++)
			{

				if(solid[y*size_x+x]==0) 
				{
					ux[y*size_x+x]=0.;
					uy[y*size_x+x]=0.;
					rho[y*size_x+x]=0.;
					for(int n=0; n<n_connect; n++)
					{
						ux[y*size_x+x] = ux[y*size_x+x]+densities[n*size_x*size_y+y*size_x+x]*e_vec_x[n];
						uy[y*size_x+x] = uy[y*size_x+x]+densities[n*size_x*size_y+y*size_x+x]*e_vec_y[n];
						rho[y*size_x+x] = rho[y*size_x+x]+densities[n*size_x*size_y+y*size_x+x];
					}
					ux[y*size_x+x] = ux[y*size_x+x]/rho[y*size_x+x];
					uy[y*size_x+x] = uy[y*size_x+x]/rho[y*size_x+x];

					// now we compute the equilibrium densities for each site. 

					for(int n=0; n<n_connect; n++)
					{
						float e_dot_u = e_vec_x[n]*ux[y*size_x+x]+e_vec_y[n]*uy[y*size_x+x];
						float u_mod_sq = ux[y*size_x+x]*ux[y*size_x+x]+uy[y*size_x+x]*uy[y*size_x+x];
						eq_densities[n*size_x*size_y+y*size_x+x]=weights[n]*(1.+(3./c)*e_dot_u+(9./2.)*(e_dot_u/c)*(e_dot_u/c)-(3./.2)*u_mod_sq/(c*c));

						// now that we have the equilibrium densities we can update the distributions. 

						temp_densities[n*size_x*size_y+y*size_x+x] = densities[n*size_x*size_y+y*size_x+x] + w*(eq_densities[n*size_x*size_y+y*size_x+x]-densities[n*size_x*size_y+y*size_x+x]);
					}
				}
			}
		}
		for(int x=0; x<size_x; x++)
		{
			for(int y=0; y<size_y; y++)
			{
				if(solid[y*size_x+x]==0) 
				{
					for(int n=0; n<n_connect; n++)
					{
						// Now that we have the updated densities, we can stream the particles. 
						//
						// n'th connect indicates the velocity direction. 
						// suppose there was N particles moving in (e_x,e_y) direction at (x,y)
						// These particles will now go from (x,y) to (x+e_x,y+e_y) lattice site. 	
						// so all that happens in densities at (x+e_x,y+e_y) gets replaced by the one 
						// at (x,y) [for the nth direction].
						//
						// This is also where we will deal with the boundary conditions. 
						//
						// In this geometry let us consider the top and the bottom borders being 
						// made of bounce back boundaries and the left and right side with 
						// periodic boundary conditions.
						// 
						
						// with the following code we can identify the PBC along x
						int x_coord,y_coord;
						if(x-e_vec_x[n]==size_x)
						{
							x_coord=0;
							y_coord=y-e_vec_y[n];
						}
						else if(x-e_vec_x[n]==-1)
						{
							x_coord=size_x-1;
							y_coord=y-e_vec_y[n];
						}
						else
						{
							x_coord=x-e_vec_x[n];
							y_coord=y-e_vec_y[n];
						}

						if(solid[y_coord*size_x+x_coord]) // this means we are trying to stream particles from a solid 
										  // node, which is not possible. 
						{
							densities[n*size_x*size_y+y*size_x+x] = temp_densities[n_opposite[n]*size_x*size_y+y*size_x+x];
						}
						else
						{
							densities[n*size_x*size_y+y*size_x+x]=temp_densities[n*size_x*size_y+y_coord*size_x+x_coord];
						}
					}
				}
			}
		}
	}
	return 0;
}
