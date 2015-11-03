
#include "global_defines.cuh"

//#define obstacles[z][y][x] obstacles[index(z,y,x)]


void   read_Reactor(const int &lx, const int &ly, const int &lz, const FLOATING &r_small, const FLOATING &s, int  *obstacles ){
	//setting obstacle

	/*c
				c.....read obstacle file
				c.....in this file you can enter the x-,and y-coordinates of
				c     any obstacles,wall boundaries are also defined here by
				c     adding single obstacles.

				.... Holes in the baffle
				c
				c                    o1
				c                 4o   o3
				c                    O             |y
				c                 5o   o6          |
				c     			   o2            |
				c                                 |--------->z*/






	//.....local variables

	int  x,y,z;
	const int xbaffle=59;
//	int iobst , m , space , d;
	int zc, yc;
//	int ycount , xcount , z_space;

	int yr,  zr;
//	int yrp, zrp;
	int R1,  R_big;
	FLOATING yshift, zshift;

	int yr1,  zr1, yctp1, zctp1;
	int yr2,  zr2, yctp2, zctp2;
	int yr3,  zr3, yctp3, zctp3;
	int yr4,  zr4, yctp4, zctp4;
	int yr5,  zr5, yctp5, zctp5;
	int yr6,  zr6, yctp6, zctp6;

	//	common /wall/ obst

	//...    tubular position and size
	R_big = 35; // in Fortran 35. it remains 35 as it doesn't participate in matrix looping
	yc= (ly +1)/2-1;//CHANGED! ORIGINALLY IT WAS yc= (ly +1)/2; AND zc= (ly +1)/2;
	zc= (ly +1)/2-1;

	// gemizei olo to domain me 0
	for ( z =0 ; z < lz ; ++z){
		for (y = 0; y < ly ; ++y){
			for ( x = 0 ; x< lx ; ++x){
				obstacles[index(z,y,x)] = 0;
			}
		}
	}

	//Big pipe
	for ( z =0 ; z < lz ; ++z){
		for (y = 0; y < ly ; ++y){
			for ( x = 0 ; x< lx ; ++x){

				zr = z - zc;
				yr = y - yc;

				//R1=pow(yr,2)+pow(zr,2);
				R1=yr*yr + zr*zr;

				if(R1 >= R_big*R_big ){
					//vazei 1 eksw apo th geometria
					obstacles[index(z,y,x)] = 1;			//original

				}

			}
		}
	}


	//baffle position
	for ( z =0 ; z < lz ; ++z){
		for (y = 0; y < ly ; ++y){
			for ( x = 0 ; x< xbaffle ; ++x){ // gemizei olo to domain me 0
				zr = z - zc;
				yr = y - yc;

				//R1=pow(yr,2)+pow(zr,2);
				R1=yr*yr + zr*zr;

				if( R1 < pow(r_small+2,2) ) {
					//frazei ta oria gurw apo ton kentriko swlhna gis aktina < r+2
					//sth sunexeia ksana adeiaze to eswteriko gia akrina <2
					obstacles[index(z,y,x)]= 1;			//original

				}

				if(   R1 <  pow(r_small,2)   ){
					obstacles[index(z,y,x)]= 0; 	//original

				}



			}
		}
	}

	//gia to baffle
	for (y = 0; y < ly ; ++y){
		for ( z =0 ; z < lz ; ++z){
			zr = z - zc;
			yr = y - yc;

			R1=yr*yr + zr*zr;
			//gemizei olo to baffle me 1, kai meta 8a "anoigei" trupes me 0 sta shmeia
			if(   R1 >=  pow(r_small+1,2) ){
				///		obstacles[z][y][xbaffle] = 1;//original

				obstacles[index(z,y,xbaffle)]= 1;
			}


			///////////////////////////////////////////
			//...  hole 1
			yctp1= (int) (yc - s);
			zctp1 = zc;
			zr1 = z - zctp1;
			yr1 = y - yctp1;
			if(yr1*yr1 + zr1*zr1 < r_small*r_small){
				///	obstacles[z][y][xbaffle] = 0;//original
				obstacles[index(z,y,xbaffle)]= 0;
			}
			////////////////////////////////////////////
			//...  hole 2
			yctp2= (int) (yc + s);
			zctp2 = zc;
			zr2 = z - zctp2;
			yr2 = y - yctp2;
			if(yr2*yr2 + zr2*zr2 < r_small*r_small){
				///	obstacles[z][y][xbaffle]= 0;//original

				obstacles[index(z,y,xbaffle)]= 0;
			}

			///////////////////////////
			zshift= s*cos(30*3.1415927/180) + 1 ;
			yshift= s*sin(30*3.1415927/180) + 1 ;
			//...  hole 3
			yctp3= yc + (int)yshift;
			zctp3 = zc - (int)zshift;
			zr3 = z - zctp3;
			yr3 = y - yctp3;
			if(yr3*yr3 + zr3*zr3 < r_small*r_small){
				///	obstacles[z][y][xbaffle] = 0;//original
				obstacles[index(z,y,xbaffle)]= 0;
			}
			//////////////////////////////////////
			//...  hole 4
			yctp4= yc + (int)yshift;
			zctp4 = zc + (int)zshift;
			zr4 = z - zctp4;
			yr4 = y - yctp4;
			if(yr4*yr4 + zr4*zr4 < r_small*r_small){
				///obstacles[z][y][xbaffle] = 0;//original
				obstacles[index(z,y,xbaffle)]= 0;

			}
			/////////////////////////////////////////
			//...  hole 5
			yctp5= yc - (int)yshift;
			zctp5 = zc + (int)zshift;
			zr5 = z - zctp5;
			yr5 = y - yctp5;
			if(yr5*yr5 + zr5*zr5 < r_small*r_small){
				///obstacles[z][y][xbaffle] = 0;//original
				obstacles[index(z,y,xbaffle)]= 0;

			}
			////////////////////////////////////////////
			//...  hole 6
			yctp6= yc-(int)yshift;
			zctp6 = zc - (int)zshift;
			zr6 = z - zctp6;
			yr6 = y - yctp6;
			if(yr6*yr6 + zr6*zr6 < r_small*r_small){
				///obstacles[z][y][xbaffle] = 0;//original
				obstacles[index(z,y,xbaffle)]= 0;

			}
		}
	}




	int error_undefined=0;
	for ( z =0 ; z < lz ; ++z){
		for (y = 0; y < ly ; ++y){
			for ( x = 0 ; x< lx ; ++x){
				if( obstacles[index(z,y,x)]!=obstacles[index(z,y,x)]){
					++error_undefined;
					cout << " BAD assignment @" << x << " " << y  << " " <<z  << "==" << error_undefined << endl;
				}
			}
		}
	}





	cout << "total errors(undefined)=" << error_undefined << endl;


}

void LBM::read_reactor(){
	/*c
					c.....read obstacle file
					c.....in this file you can enter the x-,and y-coordinates of
					c     any obstacles,wall boundaries are also defined here by
					c     adding single obstacles.

					.... Holes in the baffle
					c
					c                    o1
					c                 4o   o3
					c                    O             |y
					c                 5o   o6          |
					c     			   o2            |
					c                                 |--------->z*/
	int  x,y,z;
	const int xbaffle=59;
//	int iobst , m , space , d;
	int zc, yc;
//	int ycount , xcount , z_space;

	int yr,  zr;
//	int yrp, zrp;
	int R1,  R_big;
	FLOATING yshift, zshift;

	int yr1,  zr1, yctp1, zctp1;
	int yr2,  zr2, yctp2, zctp2;
	int yr3,  zr3, yctp3, zctp3;
	int yr4,  zr4, yctp4, zctp4;
	int yr5,  zr5, yctp5, zctp5;
	int yr6,  zr6, yctp6, zctp6;

	//	common /wall/ obst

	//...    tubular position and size
	R_big = 35; // in Fortran 35. it remains 35 as it doesn't participate in matrix looping
	yc= (ly +1)/2-1;//CHANGED! ORIGINALLY IT WAS yc= (ly +1)/2; AND zc= (ly +1)/2;
	zc= (ly +1)/2-1;

	// gemizei olo to domain me 0
	for ( z =0 ; z < lz ; ++z){
		for (y = 0; y < ly ; ++y){
			for ( x = 0 ; x< lx ; ++x){
				obstacles[index(z,y,x)] = 0;
			}
		}
	}

	//Big pipe
	for ( z =0 ; z < lz ; ++z){
		for (y = 0; y < ly ; ++y){
			for ( x = 0 ; x< lx ; ++x){

				zr = z - zc;
				yr = y - yc;

				//R1=pow(yr,2)+pow(zr,2);
				R1=yr*yr + zr*zr;

				if(R1 >= R_big*R_big ){
					//vazei 1 eksw apo th geometria
					obstacles[index(z,y,x)] = 1;			//original

				}

			}
		}
	}


	//baffle position
	for ( z =0 ; z < lz ; ++z){
		for (y = 0; y < ly ; ++y){
			for ( x = 0 ; x< xbaffle ; ++x){ // gemizei olo to domain me 0
				zr = z - zc;
				yr = y - yc;

				//R1=pow(yr,2)+pow(zr,2);
				R1=yr*yr + zr*zr;

				if( R1 < pow(r_small+2,2) ) {
					//frazei ta oria gurw apo ton kentriko swlhna gis aktina < r+2
					//sth sunexeia ksana adeiaze to eswteriko gia akrina <2
					obstacles[index(z,y,x)]= 1;			//original

				}

				if(   R1 <  pow(r_small,2)   ){
					obstacles[index(z,y,x)]= 0; 	//original

				}



			}
		}
	}

	//gia to baffle
	for (y = 0; y < ly ; ++y){
		for ( z =0 ; z < lz ; ++z){
			zr = z - zc;
			yr = y - yc;

			R1=yr*yr + zr*zr;
			//gemizei olo to baffle me 1, kai meta 8a "anoigei" trupes me 0 sta shmeia
			if(   R1 >=  pow(r_small+1,2) ){
				///		obstacles[z][y][xbaffle] = 1;//original

				obstacles[index(z,y,xbaffle)]= 1;
			}


			///////////////////////////////////////////
			//...  hole 1
			yctp1= (int) (yc - s);
			zctp1 = zc;
			zr1 = z - zctp1;
			yr1 = y - yctp1;
			if(yr1*yr1 + zr1*zr1 < r_small*r_small){
				///	obstacles[z][y][xbaffle] = 0;//original
				obstacles[index(z,y,xbaffle)]= 0;
			}
			////////////////////////////////////////////
			//...  hole 2
			yctp2= (int) (yc + s);
			zctp2 = zc;
			zr2 = z - zctp2;
			yr2 = y - yctp2;
			if(yr2*yr2 + zr2*zr2 < r_small*r_small){
				///	obstacles[z][y][xbaffle]= 0;//original

				obstacles[index(z,y,xbaffle)]= 0;
			}

			///////////////////////////
			zshift= s*cos(30*3.1415927/180) + 1 ;
			yshift= s*sin(30*3.1415927/180) + 1 ;
			//...  hole 3
			yctp3= yc + (int)yshift;
			zctp3 = zc - (int)zshift;
			zr3 = z - zctp3;
			yr3 = y - yctp3;
			if(yr3*yr3 + zr3*zr3 < r_small*r_small){
				///	obstacles[z][y][xbaffle] = 0;//original
				obstacles[index(z,y,xbaffle)]= 0;
			}
			//////////////////////////////////////
			//...  hole 4
			yctp4= yc + (int)yshift;
			zctp4 = zc + (int)zshift;
			zr4 = z - zctp4;
			yr4 = y - yctp4;
			if(yr4*yr4 + zr4*zr4 < r_small*r_small){
				///obstacles[z][y][xbaffle] = 0;//original
				obstacles[index(z,y,xbaffle)]= 0;

			}
			/////////////////////////////////////////
			//...  hole 5
			yctp5= yc - (int)yshift;
			zctp5 = zc + (int)zshift;
			zr5 = z - zctp5;
			yr5 = y - yctp5;
			if(yr5*yr5 + zr5*zr5 < r_small*r_small){
				///obstacles[z][y][xbaffle] = 0;//original
				obstacles[index(z,y,xbaffle)]= 0;

			}
			////////////////////////////////////////////
			//...  hole 6
			yctp6= yc-(int)yshift;
			zctp6 = zc - (int)zshift;
			zr6 = z - zctp6;
			yr6 = y - yctp6;
			if(yr6*yr6 + zr6*zr6 < r_small*r_small){
				///obstacles[z][y][xbaffle] = 0;//original
				obstacles[index(z,y,xbaffle)]= 0;

			}
		}
	}

}


void LBM::create_reactor_geometry_in_the_host(){
	/*c
	c.....read obstacle file
	c.....in this file you can enter the x-,and y-coordinates of
	c     any obstacles,wall boundaries are also defined here by
	c     adding single obstacles.

	.... Holes in the baffle
	c
	c                    o1
	c                 4o   o3
	c                    O             |y
	c                 5o   o6          |
	c     			   o2            |
	c                                 |--------->z*/
	int  x,y,z;
	const int xbaffle=baffle;
//	int iobst , m , space , d;
	int zc, yc;
//	int ycount , xcount , z_space;

	int yr,  zr;
//	int yrp, zrp;
	int R1,  R_big;
	FLOATING yshift, zshift;

	int yr1,  zr1, yctp1, zctp1;
	int yr2,  zr2, yctp2, zctp2;
	int yr3,  zr3, yctp3, zctp3;
	int yr4,  zr4, yctp4, zctp4;
	int yr5,  zr5, yctp5, zctp5;
	int yr6,  zr6, yctp6, zctp6;

	//	common /wall/ obst

	//...    tubular position and size
	R_big = 35; // in Fortran 35. it remains 35 as it doesn't participate in matrix looping
	yc= (ly +1)/2;//CHANGED! ORIGINALLY IT WAS yc= (ly +1)/2; AND zc= (ly +1)/2;
	zc= (lz +1)/2;


	//Big pipe
	for ( z =0 ; z < lz ; ++z){
		for (y = 0; y < ly ; ++y){
			for ( x = 0 ; x< lx ; ++x){

				zr = z - zc+1;
				yr = y - yc+1;

				//R1=pow(yr,2)+pow(zr,2);
				R1=yr*yr + zr*zr;

				if(R1 >= R_big*R_big ){
					//vazei 1 eksw apo th geometria
					obstacles[index(z,y,x)] = 1;			//original
				}
			}
		}
	}


	//baffle position
	for ( z =0 ; z < lz ; ++z){
		for (y = 0; y < ly ; ++y){
			for ( x = 0 ; x< xbaffle ; ++x){ // gemizei olo to domain me 0
				zr = z - zc+1;
				yr = y - yc+1;

				//R1=pow(yr,2)+pow(zr,2);
				R1=yr*yr + zr*zr;

				if( R1 < pow(r_small+2,2) ) {
					//frazei ta oria gurw apo ton kentriko swlhna gis aktina < r+2
					//sth sunexeia ksana adeiaze to eswteriko gia akrina <2
					obstacles[index(z,y,x)]= 1;			//original
				}

				if(   R1 <  pow(r_small,2)   ){
					obstacles[index(z,y,x)]= 0; 	//original
				}
			}
		}
	}

	//gia to baffle
	for (y = 0; y < ly ; ++y){
		for ( z =0 ; z < lz ; ++z){
			zr = z - zc+1;
			yr = y - yc+1;

			R1=yr*yr + zr*zr;
			//gemizei olo to baffle me 1, kai meta 8a "anoigei" trupes me 0 sta shmeia
			if(   R1 >=  pow(r_small+1,2) ){
				///		obstacles[z][y][xbaffle] = 1;//original
				obstacles[index(z,y,xbaffle)]= 1;
			}
			///////////////////////////////////////////
			//...  hole 1
			yctp1= (int)( yc - s);
			zctp1 = zc;
			zr1 = z - zctp1+1;
			yr1 = y - yctp1+1;
			if(yr1*yr1 + zr1*zr1 < r_small*r_small){
				///	obstacles[z][y][xbaffle] = 0;//original
				obstacles[index(z,y,xbaffle)]= 0;
			}
			////////////////////////////////////////////
			//...  hole 2
			yctp2= (int)( yc + s);
			zctp2 = zc;
			zr2 = z - zctp2+1;
			yr2 = y - yctp2+1;
			if(yr2*yr2 + zr2*zr2 < r_small*r_small){
				///	obstacles[z][y][xbaffle]= 0;//original
				obstacles[index(z,y,xbaffle)]= 0;
			}
			///////////////////////////
			zshift= s*cos(30*3.1415927/180) + 1 ;
			yshift= s*sin(30*3.1415927/180) + 1 ;
			//...  hole 3
			yctp3= yc + (int)(yshift);
			zctp3 = zc - (int)(zshift);
			zr3 = z - zctp3+1;
			yr3 = y - yctp3+1;
			if(yr3*yr3 + zr3*zr3 < r_small*r_small){
				///	obstacles[z][y][xbaffle] = 0;//original
				obstacles[index(z,y,xbaffle)]= 0;
			}
			//////////////////////////////////////
			//...  hole 4
			yctp4= yc + (int)(yshift);
			zctp4 = zc + (int)(zshift);
			zr4 = z - zctp4+1;
			yr4 = y - yctp4+1;
			if(yr4*yr4 + zr4*zr4 < r_small*r_small){
				///obstacles[z][y][xbaffle] = 0;//original
				obstacles[index(z,y,xbaffle)]= 0;
			}
			/////////////////////////////////////////
			//...  hole 5
			yctp5= yc - (int)(yshift);
			zctp5 = zc + (int)(zshift);
			zr5 = z - zctp5+1;
			yr5 = y - yctp5+1;
			if(yr5*yr5 + zr5*zr5 < r_small*r_small){
				///obstacles[z][y][xbaffle] = 0;//original
				obstacles[index(z,y,xbaffle)]= 0;
			}
			////////////////////////////////////////////
			//...  hole 6
			yctp6= yc-(int)(yshift);
			zctp6 = zc - (int)(zshift);
			zr6 = z - zctp6+1;
			yr6 = y - yctp6+1;
			if(yr6*yr6 + zr6*zr6 < r_small*r_small){
				///obstacles[z][y][xbaffle] = 0;//original
				obstacles[index(z,y,xbaffle)]= 0;
			}
		}
	}
}
