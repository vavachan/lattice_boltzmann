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
#include<fstream>

int main()
{
	// these are the sizes along x direction and y direction, say in meters
	const float lengthX=2.2;
	const float lengthY=.41;

	const float flow_vel = 0.5; // in m/s so it takes about 11 seconds to flow from one end 
				    // the other.
        const float charphysL = 0.05; // this is the radius of 
					// the cylinder. its a tiny tiny cylinder. 
					//
	const float viscosity=0.2e-3; // kinematic viscosity m^2/s;
				    // I think this is the viscosity of water. 
				    //
	std::cout<<charphysL*flow_vel/viscosity<<" = Re \n";
	// The above values represent what is given to us 
	// as the system details 
	//
	// The reynolds number can be calculated. 
	//
	// 	Re = charL * flow_vel / viscosity \approx 100
	//
	//
	//const float Re = 10.; // Reynolds number, how do we fix for diffusion ? I don't know yet
        const float N = 20.; // Resolution - this is upto us 
	const float tau = 0.53; // Relaxation time - this is also upto us 
				//
	// The above two variables needs to changed to get a stable simulation. 

	const float L = charphysL/N; // this is also delta_x
	const float c =1.;// sqrt(1./3.); // this is speed of sound ??!! in lattice units.
			  // we will worry about this later. 

	const float delta_t = 1./3.*(tau-0.5)*L*L/viscosity; 

	std::cout<<flow_vel/(L/delta_t)<<"\n";
	// average velocity is towards the right 
	
	float uox=flow_vel/(L/delta_t);
	float uoy=0.;

	float cylinder_radius = charphysL/L;
	int timesteps = 1000;
	int write_after = 100;

	std::cout<<L/delta_t<<"\t"<<delta_t<<"\t"<<cylinder_radius<<"\n";;
	
	const int size_x = (lengthX+L)/L;
	const int size_y = (lengthY+L)/L; 

	std::cout<<size_x<<"\t"<<size_y<<"\t"<<size_x*size_y<<"\n";

	//float m=1.;
	//float k=1;	
	//float T=1.;
			       //
	//std::cout<<c<<"\n";

	float w = 1./tau; // this is related to updating of the distribution in the collision step

	// lattice connectivity 
	int n_connect = 9;

	float * densities = nullptr;
	
	densities = new float[n_connect*size_x*size_y];

	float * eq_densities = nullptr;
	
	eq_densities = new float[n_connect*size_x*size_y];

	float * temp_densities = nullptr;
	temp_densities = new float[n_connect*size_x*size_y]; // I need these. 
							  
	int * solid = nullptr;

	solid = new int[size_x*size_y];
	// initialization
	
	for(int x=0; x<size_x; x++)
	{
		for(int y=0; y<size_y; y++)
		{
		//	if(y==size_y-1 or y==0)
		//	{
		//		solid[y*size_x+x]=1.;
		//	}
		//	else 
				solid[y*size_x+x]=0.;
		}
	}
	
	// insert a circular disk as a obstacle. 
	
	// we will position the disk at x=5, y=4 
	// with a radius 1. 
	int xo=int(1./5.*size_x);
	int yo=int(1./2.*size_y);
	int n_solid=0;
	for(int x=0; x<size_x; x++)
	{
		for(int y=0; y<size_y; y++)
		{
			if((x-xo)*(x-xo)+(y-yo)*(y-yo)<=pow(cylinder_radius,2))
			{
				//std::cout<<x<<"\t"<<y<<"\n";
				solid[y*size_x+x]=1.;
				n_solid = n_solid+1;

			}	
		}
	}
	
//	for(int y=0;y<size_y;y++)
//	{
//		solid[y*size_x+size_x/2]=1.;
//	}

	std::cout<<n_solid<<"\n";
	// n_connect is the number of velocities in the system. 
	
	// the following are the vectors which corresponds to the 
	// n_connect velocities. 

	int * e_vec_x = new int[n_connect];
	int * e_vec_y = new int[n_connect];
	
	int * n_opposite = new int[n_connect];

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


	for(int x=0; x<size_x; x++)
	{
		//int x=0;
		for(int y=0; y<size_y; y++)
		{
			if(solid[y*size_x+x]==0)
			{
				for(int n=0; n<n_connect; n++)
				{
					float e_dot_u = e_vec_x[n]*uox+e_vec_y[n]*uoy;
					float u_mod_sq = uox*uox+uoy*uoy;
					//if(x==size_x/2)
					{
						densities[n*size_x*size_y+y*size_x+x]=weights[n]*(1.+3.*e_dot_u+(9./2.)*(e_dot_u)*(e_dot_u)-(3./2.)*u_mod_sq);
					}
				//	else
				//	{
				//		densities[n*size_x*size_y+y*size_x+x]=weights[n];
				//	}

					//std::cout<<n<<"\t"<<weights[n]*(1.+(3./c)*e_dot_u+(9./2.)*(e_dot_u/c)*(e_dot_u/c)-(3./2.)*u_mod_sq/(c*c))<<"\n";
				}
			}
		}
	}
	//
	// we need a few more arrays to store the average velocity 
	// and the density

	float * ux=new float[size_x*size_y];
	float * uy=new float[size_x*size_y];
	float * rho=new float[size_x*size_y];
	// now we do the simulation 


	std::fstream a_vel;
	std::fstream den_file;
	std::string filename;

	//std::cout<<t<<"\n";
	int t=-1;
	filename = "average_velocity_"+std::to_string(t)+".dat";
	a_vel.open(filename,std::ios::out);
	filename = "densities_"+std::to_string(t)+".dat";
	den_file.open(filename,std::ios::out);
	for(int x=0; x<size_x; x++)
	{
		for(int y=0; y<size_y; y++)
		{
			den_file<<x<<"\t"<<y<<"\t";
			for(int n=0;n<n_connect; n++)
			{
				den_file<<densities[n*size_x*size_y+y*size_x+x]-weights[n]<<"\t";
			}
			den_file<<"\n";
			a_vel<<x<<"\t"<<y<<"\t"<<sqrt(ux[y*size_x+x]*ux[y*size_x+x]+uy[y*size_x+x]*uy[y*size_x+x])<<"\t"<<ux[y*size_x+x]<<"\t"<<uy[y*size_x+x]<<"\n";
		}
	}
	a_vel.close();
	den_file.close();
	int domain=size_x*size_y;
	for(int t=0; t<timesteps;t++)
	{
		for(int x=0; x<size_x; x++)
		{
			for(int y=0; y<size_y; y++)
			{

				int pos=y*size_x+x;
				ux[pos]=0.;
				uy[pos]=0.;
				rho[pos]=0.;
				if(solid[pos]==0) 
				{

					for(int n=0; n<n_connect; n++)
					{
						ux[pos] = ux[pos]+densities[n*domain+pos]*e_vec_x[n];
						uy[pos] = uy[pos]+densities[n*domain+pos]*e_vec_y[n];
						rho[pos] = rho[pos]+densities[n*domain+pos];
					}
					ux[pos] = ux[pos]/rho[pos];
					uy[pos] = uy[pos]/rho[pos];


					// now we compute the equilibrium densities for each site. 

					for(int n=0; n<n_connect; n++)
					{
						float e_dot_u = e_vec_x[n]*ux[pos]+e_vec_y[n]*uy[pos];
						float u_mod_sq = ux[pos]*ux[pos]+uy[pos]*uy[pos];

						// now that we have the equilibrium densities we can update the distributions. 

						temp_densities[n*domain+pos] = densities[n*domain+pos] + w*(rho[pos]*weights[n]*(1.+3.*e_dot_u+(9./2.)*(e_dot_u)*(e_dot_u)-(3./2.)*u_mod_sq)-densities[n*domain+pos]);
					}
				}
				else
				{
					// this is a solid node.
					// here for each n we replace it with the opposite direction.
					for(int n=0; n<n_connect; n++)
					{
						//temp_densities[n_opposite[n]*domain+pos] = densities[n*domain+pos];
						temp_densities[n*domain+pos] = densities[n_opposite[n]*domain+pos];
					}
				}
			}
		}
	  //	for(int n=0; n<n_connect; n++)
	  //	{
	  //		//for(int x=0; x<size_x; x++)
	  //		{
	  //			for(int y=0; y<size_y; y++)
	  //			{
	  //				float e_dot_u = e_vec_x[n]*uox+e_vec_y[n]*uoy;
	  //				float u_mod_sq = uox*uox+uoy*uoy;
	  //				temp_densities[n*domain+y*size_x]=weights[n]*(1.+(3./c)*e_dot_u+(9./2.)*(e_dot_u/c)*(e_dot_u/c)-(3./2.)*u_mod_sq/(c*c));
	  //			}
	  //		}
	  //	}

		// This is the propogation step
		for(int x=0; x<size_x; x++)
		{
			for(int y=0; y<size_y; y++)
			{
				//if(solid[y*size_x+x]==0) 
				int pos=y*size_x+x;
				{
					for(int n=0; n<n_connect; n++)
					{
						// 		* This was written for half-bounce-back**
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

						/*
						 * 		THIS IS FOR Full_BOUNCE_BACK
						 *
						 * Here the position we are in right now - (x,y) can be solid 
						 * or not. 
						 *
						 * Even a solid node contains populations. 
						 *
						 * So we will update the densities of the neighbouring nodes here. 
						 *
						 * i.e. We will update (x+en_x,y+en_y) densities here. 
						 */
						
						// with the following code we can identify the PBC along x
						//int x_coord,y_coord; // x_coord and y_coord represent the coordinates from which we will stream the particles to x,y

						int x_coord,y_coord; // x_coord and y_coord represent the coordinates **to** which we will stream the particles **from** x,y
						if(x+e_vec_x[n]==size_x)
						{
							x_coord=0;
							//y_coord=y-e_vec_y[n];
						}
						else if(x+e_vec_x[n]==-1)
						{
							x_coord=size_x-1;
							//y_coord=y-e_vec_y[n];
						}
						else
						{
							x_coord=x+e_vec_x[n];
							//y_coord=y-e_vec_y[n];
						}
						if(y+e_vec_y[n]==size_y)
						{
							y_coord=0;
							//y_coord=y-e_vec_y[n];
						}
						else if(y+e_vec_y[n]==-1)
						{
							y_coord=size_y-1;
							//y_coord=y-e_vec_y[n];
						}
						else
						{
							y_coord=y+e_vec_y[n];
							//y_coord=y-e_vec_y[n];
						}

						//if(solid[y_coord*size_x+x_coord]) // this means we are trying to stream particles from a solid 
						//				  // node, which is not possible. 
						//{
						//	densities[n*domain+y*size_x+x] = temp_densities[n_opposite[n]*domain+y*size_x+x]; // here we are implementing the half-wall reflection 
						//													// boundary condition 
						//													// Since there are not particles that are streaming from 
						//													// the solid wall, we don't have a value for densities 
						//													// for n_connects from the wall to x,y. 
						//													// So we put those densitiy values as particles which got 
						//													// reflected off the wall (...)
						//}
						//else
					//	{
					//		//densities[n*domain+y*size_x+x]=temp_densities[n*domain+y_coord*size_x+x_coord];
					//	}
						densities[n*domain+y_coord*size_x+x_coord]=temp_densities[n*domain+pos];
					}
				}
			}
		}
		if(t%write_after==0)
		{
			std::cout<<t<<"\n";
			filename = "average_velocity_"+std::to_string(t)+".dat";
			a_vel.open(filename,std::ios::out);
			filename = "densities_"+std::to_string(t)+".dat";
			den_file.open(filename,std::ios::out);
			for(int x=0; x<size_x; x++)
			{
				for(int y=0; y<size_y; y++)
				{
					den_file<<x*(L)<<"\t"<<y*L<<"\t";
					for(int n=0;n<n_connect; n++)
					{
						den_file<<densities[n*domain+y*size_x+x]-weights[n]<<"\t";
					}
					den_file<<"\n";
					a_vel<<x<<"\t"<<y<<"\t"<<sqrt(ux[y*size_x+x]*ux[y*size_x+x]+uy[y*size_x+x]*uy[y*size_x+x])<<"\t"<<ux[y*size_x+x]<<"\t"<<uy[y*size_x+x]<<"\n";
				}
			}
			a_vel.close();
			den_file.close();
		}
	}
	//for(int x=0; x<size_x; x++)
	//{
	//	for(int y=0; y<size_y; y++)
	//	{
	//		std::cout<<x<<"\t"<<y<<"\t"<<densities[y*size_x+x]<<"\t"<<solid[y*size_x+x]<<"\n";
	//	}
	//}
	return 0;
}
