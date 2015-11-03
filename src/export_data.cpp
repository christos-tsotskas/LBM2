#include "global_defines.cuh"

#define STR1(x)  #x
#define STRINGIFY_MACRO(x)  STR1(x)

//void  LBM::save_data(   FLOATING &pr_diff, FLOATING &vor, const int& iteration,
//		FLOATING *Ux, FLOATING *Uy, FLOATING *Uz, FLOATING *Pressure, FLOATING *Wx, FLOATING *Wy, FLOATING *Wz, const FLOATING &time_elapsed){
FLOATING LBM::calculate_vorticity(){
	FLOATING vor=0.0;
	FLOATING w_x,w_y,w_z,w_out ;


	for (int z = 0 ; z< lz ; ++z){
		for (int y = 0 ; y< ly ; ++y){
			for (int x = 0 ; x< lx ; ++x){
				if (obstacles[index(z,y,x)]==1) {
					w_x = 0.0;
					w_y = 0.0;
					w_z = 0.0;

				}else if(x==0){
					w_x=0.125*((Uy[index(z-1,y-1,x)]+2*Uy[index(z-1,y,x)]+Uy[index(z-1,y+1,x)])+
							(Uz[index(z-1,y+1,x)]+2*Uz[index(z,y+1,x)]+Uz[index(z+1,y+1,x)])-
							(Uy[index(z+1,y+1,x)]+2*Uy[index(z+1,y,x)]+Uy[index(z+1,y-1,x)])-
							(Uz[index(z+1,y-1,x)]+2*Uz[index(z,y-1,x)]+Uz[index(z-1,y-1,x)]));

					w_y=0.125*((Uz[index(z-1,y,x)]+2*Uz[index(z,y,x)]+Uz[index(z+1,y,x)])+
							(Ux[index(z+1,y,x)]+2*Ux[index(z+1,y,x+1)]+Ux[index(z+1,y,x+2)])-
							(Uz[index(z+1,y,x+2)]+2*Uz[index(z,y,x+2)]+Uz[index(z-1,y,x+2)])-
							(Ux[index(z-1,y,x+2)]+2*Ux[index(z-1,y,x+1)]+Ux[index(z-1,y,x)]));

					w_z=0.125*((Ux[index(z,y-1,x)]+2*Ux[index(z,y-1,x+1)]+Ux[index(z,y-1,x+2)])+
							(Uy[index(z,y-1,x+2)]+2*Uy[index(z,y,x+2)]+Uy[index(z,y+1,x+2)])-
							(Ux[index(z,y+1,x+2)]+2*Ux[index(z,y+1,x+1)]+Ux[index(z,y+1,x)])-
							(Uy[index(z,y+1,x)]+2*Uy[index(z,y,x)]+Uy[index(z,y-1,x)]));

				}else if(x==(lx-1)){
					w_x=0.125*((Uy[index(z-1,y-1,x)]+2*Uy[index(z-1,y,x)]+Uy[index(z-1,y+1,x)])+
							(Uz[index(z-1,y+1,x)]+2*Uz[index(z,y+1,x)]+Uz[index(z+1,y+1,x)])-
							(Uy[index(z+1,y+1,x)]+2*Uy[index(z+1,y,x)]+Uy[index(z+1,y-1,x)])-
							(Uz[index(z+1,y-1,x)]+2*Uz[index(z,y-1,x)]+Uz[index(z-1,y-1,x)]));

					w_y=0.125*((Uz[index(z-1,y,x-2)]+2*Uz[index(z,y,x-2)]+Uz[index(z+1,y,x-2)])+
							(Ux[index(z+1,y,x-2)]+2*Ux[index(z+1,y,x-1)]+Ux[index(z+1,y,x)])-
							(Uz[index(z+1,y,x)]+2*Uz[index(z,y,x)]+Uz[index(z-1,y,x)])-
							(Ux[index(z-1,y,x)]+2*Ux[index(z-1,y,x-1)]+Ux[index(z-1,y,x-2)]));

					w_z=0.125*((Ux[index(z,y-1,x-2)]+2*Ux[index(z,y-1,x-1)]+Ux[index(z,y-1,x)])+
							(Uy[index(z,y-1,x)]+2*Uy[index(z,y,x)]+Uy[index(z,y+1,x)])-
							(Ux[index(z,y+1,x)]+2*Ux[index(z,y+1,x-1)]+Ux[index(z,y+1,x-2)])-
							(Uy[index(z,y+1,x-2)]+2*Uy[index(z,y,x-2)]+Uy[index(z,y-1,x-2)]));
				}else{
					w_x=0.125*((Uy[index(z-1,y-1,x)]+2*Uy[index(z-1,y,x)]+Uy[index(z-1,y+1,x)])+
							(Uz[index(z-1,y+1,x)]+2*Uz[index(z,y+1,x)]+Uz[index(z+1,y+1,x)])-
							(Uy[index(z+1,y+1,x)]+2*Uy[index(z+1,y,x)]+Uy[index(z+1,y-1,x)])-
							(Uz[index(z+1,y-1,x)]+2*Uz[index(z,y-1,x)]+Uz[index(z-1,y-1,x)]));

					w_y=0.125*((Uz[index(z-1,y,x-1)]+2*Uz[index(z,y,x-1)]+Uz[index(z+1,y,x-1)])+
							(Ux[index(z+1,y,x-1)]+2*Ux[index(z+1,y,x)]+Ux[index(z+1,y,x+1)])-
							(Uz[index(z+1,y,x+1)]+2*Uz[index(z,y,x+1)]+Uz[index(z-1,y,x+1)])-
							(Ux[index(z-1,y,x+1)]+2*Ux[index(z-1,y,x)]+Ux[index(z-1,y,x-1)]));

					w_z=0.125*((Ux[index(z,y-1,x-1)]+2*Ux[index(z,y-1,x)]+Ux[index(z,y-1,x+1)])+
							(Uy[index(z,y-1,x+1)]+2*Uy[index(z,y,x+1)]+Uy[index(z,y+1,x+1)])-
							(Ux[index(z,y+1,x+1)]+2*Ux[index(z,y+1,x)]+Ux[index(z,y+1,x-1)])-
							(Uy[index(z,y+1,x-1)]+2*Uy[index(z,y,x-1)]+Uy[index(z,y-1,x-1)]));
				}
				Wx[index(z,y,x)]=w_x;
				Wy[index(z,y,x)]=w_y;
				Wz[index(z,y,x)]=w_z;

				w_out=sqrt(w_x*w_x+w_y*w_y+w_z*w_z);
				vor+=w_out;
			}
		}
	}

	return vor;
}

FLOATING LBM::calculate_pressure_differences(){
	FLOATING  rho=0.0,press=0.0;
	FLOATING  pr_out=0.0,pr_in=0.0;

	for (int z = 0 ; z< lz ; ++z){
		for (int y = 0 ; y< ly ; ++y){
			for (int x = 0 ; x< lx ; ++x){


				//.........if obstacle node, nothing is to do ...
				if (  obstacles[index(z,y,x)]==1 ) {
					//...........obstacle indicator
					//...........velocity components = 0
					Ux[index(z,y,x)] = 0.0;
					Uy[index(z,y,x)] = 0.0;
					Uz[index(z,y,x)] = 0.0;
					//..........Temperature
					//!T_t = T[z][y][x]
					//...........pressure = average pressure
					press = density * c_squ;


				}else{
					//...........integral local density
					//...........initialize variable ro
					rho=0.0;
					rho+=D3.Q0[index(z,y,x)]+D3.Q1[index(z,y,x)]+D3.Q2[index(z,y,x)]+D3.Q3[index(z,y,x)];
					rho+=D3.Q4[index(z,y,x)]+D3.Q5[index(z,y,x)]+D3.Q6[index(z,y,x)]+D3.Q7[index(z,y,x)];
					rho+=D3.Q8[index(z,y,x)]+D3.Q9[index(z,y,x)]+D3.Q10[index(z,y,x)]+D3.Q11[index(z,y,x)];
					rho+=D3.Q12[index(z,y,x)]+D3.Q13[index(z,y,x)]+D3.Q14[index(z,y,x)]+D3.Q15[index(z,y,x)];
					rho+=D3.Q16[index(z,y,x)]+D3.Q17[index(z,y,x)]+D3.Q18[index(z,y,x)];



					//...........x-, and y- velocity components
					//FOLLOW DJENIDI NOTATION!!!!


					Ux[index(z,y,x)] = (FLOATING) (D3.Q1[index(z,y,x)] + D3.Q7[index(z,y,x)] + D3.Q10[index(z,y,x)] +
							D3.Q11[index(z,y,x)] + D3.Q14[index(z,y,x)] -
							(D3.Q3[index(z,y,x)] + D3.Q8[index(z,y,x)] + D3.Q9[index(z,y,x)] +
									D3.Q12[index(z,y,x)] + D3.Q13[index(z,y,x)])) / rho;

					Uy[index(z,y,x)] = (FLOATING) (D3.Q2[index(z,y,x)] + D3.Q8[index(z,y,x)] + D3.Q7[index(z,y,x)] +
							D3.Q16[index(z,y,x)] + D3.Q15[index(z,y,x)] -
							(D3.Q4[index(z,y,x)] + D3.Q9[index(z,y,x)] + D3.Q10[index(z,y,x)] +
									D3.Q17[index(z,y,x)] + D3.Q18[index(z,y,x)])) / rho;

					Uz[index(z,y,x)] = (FLOATING) (D3.Q5[index(z,y,x)] + D3.Q13[index(z,y,x)] + D3.Q14[index(z,y,x)] +
							D3.Q15[index(z,y,x)] + D3.Q18[index(z,y,x)] -
							(D3.Q6[index(z,y,x)] + D3.Q12[index(z,y,x)] + D3.Q11[index(z,y,x)] +
									D3.Q16[index(z,y,x)] + D3.Q17[index(z,y,x)])) / rho;




					//...........pressure
					press = rho * c_squ;


				}
				Pressure[index(z,y,x)]=press;


				if( x==lx-2){//penultimate slice!!!
					pr_out+=press;
					//					temp_output <<  "[" <<x << ", " << y << ", " << z << "]." <<" pr_+ " << press  << "rho ("<<rho  <<  ") * c_squ(" << c_squ << ") "<< endl;
					//					temp_output << " pr_out " << pr_out << endl;
				}

				if( x==1){//second slice!!!!!
					pr_in+=press;
					//					temp_output <<  "[" <<x << ", " << y << ", " << z << "]." << " pr_- " << press << "rho ("<<rho  <<  ") * c_squ(" << c_squ << ") "<< endl;
					//					temp_output << " pr_in " << pr_in << endl;
				}
			}
		}
	}
	return pr_in-pr_out;//// IN!
}

void LBM::write_log_file(const int iteration){
	time (&time_end);
	time_elapsed = difftime (time_end,time_start);

	char buffer1[256];

	//output file
	snprintf(buffer1, sizeof(buffer1), "LBM2_%s_report_at_iteration_%d.log",case_name.c_str(), iteration);

	ofstream log_file(buffer1);


	log_file << "Case: "<< case_name.c_str() << endl;

	log_file << "LBM Geometry:" << endl;
	log_file << "\tX:" << lx << endl;
	log_file << "\tY:" << ly << endl;
	log_file << "\tZ:" << lz << endl;

	log_file << "Configuration Parameters:" << endl;
	log_file << "\titerations: " << max_iterations << endl;
	log_file << "\tcheck step: " << check_step << endl;
	log_file << "\tnu: " << nu<< endl;
	log_file << "\tr_small: " << r_small<< endl;
	log_file << "\treynolds: " << reynolds<< endl;
	log_file << "\ts: " << s<< endl;
	log_file << "\tbaffle position on X=" << baffle << endl;
	log_file << "\tCUDA threads per kernel: " << threads_per_kernel<< endl;

	log_file << "Performance Criteria" << endl;
	log_file <<  "\t pressure difference: " << pr_diff << "(pr_in:"<< pr_in << "- pr_out:" << pr_out <<")" << endl;
	log_file <<  "\t vorticity: " << vor << endl;

	log_file << "Elapsed Time: " << time_elapsed << " seconds "<< endl;
	log_file.close();
}

void LBM::write_objectives(){
	//output file
	//objectives: for interfacing with MOTS2
	ofstream objectives( "LBM2_objectives.txt"  );
	objectives<< vor << "\t" << pr_diff << endl;
	objectives.close();
}

void  LBM::calculate_macroscopic_quantities(const int& iteration){



	//todo: check, maybe rec_pr's size might have conflict in size dimension
	/************************************************************************
	 *                                                                      *
	 *     Output of results to file 'tecplot.dat'                          *
	 *                                                                      *
	 *            Christos Tsotskas											 *
	 *                                                                      *
	 *     Last change: 09/06/2013                                          *
	 *                                                                      *
	 ************************************************************************/

	if(data_location==GPU)
		copy_data_from_device_to_host();

	pr_diff = calculate_pressure_differences();

	vor = calculate_vorticity();

	write_log_file(iteration);
	write_objectives();

	cout << " save_data exit" << endl;
}

void  LBM::write_VTK_SI(const int& iteration){

	//void  LBM::write_VTK_SI(const int& iteration,
	//		const FLOATING *Ux, const FLOATING *Uy, const FLOATING *Uz,
	//		const FLOATING *Pressure, const FLOATING *Wx, const FLOATING *Wy, const FLOATING *Wz){
	//todo: check, maybe rec_pr's size might have conflict in size dimension


	if(data_location==GPU)
		copy_data_from_device_to_host();
	//.....local variables
	int  x,y,z;
	FLOATING  pr_out=0.0,pr_in=0.0;
	FLOATING  rho=0.0,press=0.0;

	char buffer1[256];




	//.....square speed of sound

	//.....open results output file
	//.....write header for postprocessing with TECPLOT software
	/*c.....loop over all nodes
	c.....attention: actual densities are stored after the propagation
	c                step in the help-array n_hlp !
	c*/

	//output file
	snprintf(buffer1, sizeof(buffer1), "LBM2_%s_data_SI_iteration_%d.vtk", case_name.c_str(),iteration);
	ofstream VTK_file(buffer1);

	VTK_file<<"# vtk DataFile Version 3.0"<< endl;
	VTK_file<<"vtk output"<< endl;
	VTK_file<<"ASCII"<< endl;

	FLOATING conversion=12.73;
	FLOATING pressure_conversion=162400;
	int next_line_counter=0;



	//coordinates
	VTK_file<<"DATASET STRUCTURED_GRID"<< endl;
	VTK_file<<"DIMENSIONS "<< lx <<" "<< ly <<" " << lz << endl;
	VTK_file<<"POINTS " << lx*ly*lz <<" "STRINGIFY_MACRO(FLOATING)<< endl;

	for (z = 0 ; z< lz ; ++z)
		for (y = 0 ; y< ly ; ++y)
			for (x = 0 ; x< lx ; ++x){
				VTK_file << x << "\t" << y << "\t" << "\t" << z << "\t";
				++next_line_counter;
				if(next_line_counter%9==0)
					VTK_file << endl;
			}

	VTK_file<< endl << endl;

	//print Pressure
	VTK_file<<"POINT_DATA "<< lx*ly*lz << endl;
	VTK_file<<"SCALARS Pressure "STRINGIFY_MACRO(FLOATING)<< endl;
	VTK_file<<"LOOKUP_TABLE default"<< endl;
	next_line_counter=0;
	for (z = 0 ; z< lz ; ++z){
		for (y = 0 ; y< ly ; ++y){
			for (x = 0 ; x< lx ; ++x){
				//.........if obstacle node, nothing is to do ...
				VTK_file <<  Pressure[index(z,y,x)]*pressure_conversion << "\t";

				++next_line_counter;
				if(next_line_counter%9==0)
					VTK_file << endl;

			}
		}
	}
	VTK_file<< endl << endl;



	//print Velocity Vector (U,V,W)
	VTK_file<<"VECTORS Velocity "STRINGIFY_MACRO(FLOATING)<< endl;
	next_line_counter=0;
	for (z = 0 ; z< lz ; ++z){
		for (y = 0 ; y< ly ; ++y){
			for (x = 0 ; x< lx ; ++x){
				//.........if obstacle node, nothing is to do ...
				VTK_file <<  Ux[index(z,y,x)]*conversion << "\t" << Uy[index(z,y,x)]*conversion << "\t" << Uz[index(z,y,x)]*conversion << "\t";

				++next_line_counter;
				if(next_line_counter%9==0)
					VTK_file << endl;

			}
		}
	}
	VTK_file<< endl << endl;

	//print Velocity Vector (Wx,Wy,Wz)
	VTK_file<<"VECTORS Vorticity "STRINGIFY_MACRO(FLOATING)<< endl;
	next_line_counter=0;
	for (z = 0 ; z< lz ; ++z){
		for (y = 0 ; y< ly ; ++y){
			for (x = 0 ; x< lx ; ++x){
				//.........if obstacle node, nothing is to do ...
				VTK_file <<  Wx[index(z,y,x)]*conversion << "\t" << Wy[index(z,y,x)]*conversion << "\t" << Wz[index(z,y,x)]*conversion << "\t";

				++next_line_counter;
				if(next_line_counter%9==0)
					VTK_file << endl;

			}
		}
	}
	VTK_file<< endl << endl;


	//.....close file 'Tecplot.dat'
	cout << "vtk_SI exit" << endl;
	VTK_file.close();

}

void create_VTK_header(const string ){

}

void LBM::geometry_file_in_VTK( ){

	char buffer1[256];
	int x,y,z;

	snprintf(buffer1, sizeof(buffer1), "LBM2_%s_geometry.vtk", case_name.c_str());
	ofstream VTK_file(buffer1);

	VTK_file<<"# vtk DataFile Version 3.0"<< endl;
	VTK_file<<"vtk output"<< endl;
	VTK_file<<"ASCII"<< endl;

	int next_line_counter=0;

	//coordinates
	VTK_file<<"DATASET STRUCTURED_GRID"<< endl;
	VTK_file<<"DIMENSIONS "<< lx <<" "<< ly <<" " << lz << endl;
	VTK_file<<"POINTS " << lx*ly*lz <<" " STRINGIFY_MACRO(FLOATING) << endl;

	for (z = 0 ; z< lz ; ++z)
		for (y = 0 ; y< ly ; ++y)
			for (x = 0 ; x< lx ; ++x){
				VTK_file << x << "\t" << y << "\t" << "\t" << z << "\t";
				++next_line_counter;
				if(next_line_counter%9==0)
					VTK_file << endl;
			}

	VTK_file<< endl << endl;

	//print Obstacles
	VTK_file<<"POINT_DATA "<< lx*ly*lz << endl;
	VTK_file<<"SCALARS Obstacles int"<< endl;
	VTK_file<<"LOOKUP_TABLE default"<< endl;
	next_line_counter=0;
	for (z = 0 ; z< lz ; ++z){
		for (y = 0 ; y< ly ; ++y){
			for (x = 0 ; x< lx ; ++x){
				//.........if obstacle node, nothing is to do ...
				VTK_file <<  obstacles[index(z,y,x)] << "\t";

				++next_line_counter;
				if(next_line_counter%9==0)
					VTK_file << endl;
			}
		}
	}
	VTK_file<< endl << endl;

	VTK_file.close();
	cout << "geometry file exported!" << endl;

}
