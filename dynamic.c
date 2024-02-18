#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

#define WIDTH 640
#define HEIGHT 480
#define MAX_ITER 255

#define imag_min -2.0
#define imag_max 3.0
#define real_min -2.0
#define real_max 3.0

struct complex {
  double real;
  double imag;
};

//cal_pixel and save_pgm functions were not changed from the sequential code

int cal_pixel(struct complex c) {

    double z_real = 0;
    double z_imag = 0;

    double z_real2, z_imag2, lengthsq;

    int iter = 0;
    do {
        z_real2 = z_real * z_real;
        z_imag2 = z_imag * z_imag;

        z_imag = 2 * z_real * z_imag + c.imag;
        z_real = z_real2 - z_imag2 + c.real;
        lengthsq =  z_real2 + z_imag2;
        iter++;
    }
    while ((iter < MAX_ITER) && (lengthsq < 4.0));

    return iter;

}

void save_pgm(const char *filename, int image[HEIGHT][WIDTH]) {
    FILE* pgmimg; 
    int temp;
    pgmimg = fopen(filename, "wb"); 
    fprintf(pgmimg, "P2\n"); // Writing Magic Number to the File   
    fprintf(pgmimg, "%d %d\n", WIDTH, HEIGHT);  // Writing Width and Height
    fprintf(pgmimg, "255\n");  // Writing the maximum gray value 
    int count = 0; 
    
    for (int i = 0; i < HEIGHT; i++) { 
        for (int j = 0; j < WIDTH; j++) { 
            temp = image[i][j]; 
            fprintf(pgmimg, "%d ", temp); // Writing the gray values in the 2D array to the file 
        } 
        fprintf(pgmimg, "\n"); 
    } 
    fclose(pgmimg); 
} 


int main(int argc, char *argv[]) {
	//Initialize MPI Details such as processor rank and # of nodes
	int size, rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	
	//Initialize local 1D array and final_image array
	int partial_image[WIDTH];
	int final_image[HEIGHT][WIDTH];
	
	//Further initializations for dynamic allocation of work
	//Count and row initialized to zero and will be changed according to which processors are being used
	//Terminate value for when row == HEIGHT
	//Complex structure to hold real and imaginary values 
	int count, row = 0;
	int terminate = -1;
	struct complex c;
	
	//Start Timer
	double start_time = MPI_Wtime();
    	
    	//MASTER PROGRAM
    	if(rank==0){
    		//Iterate through p-1 processes and send initial row numbers.
    		//1 ROW PER PROCESSOR - increment row to keep track and count to know how many processors are currently working
    		
		for(int p = 1;p<size;p++) {
    			MPI_Send(&row,1,MPI_INT,p,0,MPI_COMM_WORLD);
    			count++;
    			row++;
  		}
  		MPI_Status status;
  		//WHILE LOOP terminates when no more processors are working
  		do {
  			//Receive 1-D row of values from any slave and decrement count
        		MPI_Recv(partial_image,WIDTH,
        			MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD, &status);
        		count--;
        		
        		//Keep track of which processor was just received and the row it sent
        		int proc = status.MPI_SOURCE;
        		int rcvdRow = status.MPI_TAG;
        		//Check if image has been completed
        		//If not, send rows to complete to other available processors
        		if(row < HEIGHT){
        			MPI_Send(&row,1,MPI_INT,proc, 0, MPI_COMM_WORLD);
        			row++;
        			count++;
        		}
        		//If image is finished, send terminate tag to let slaves break the while loop
        		else{
        			MPI_Send(&terminate,1,MPI_INT,proc,0,MPI_COMM_WORLD);
        		}
        		//Add the received row to final_image array
        		for( int col = 0;col<WIDTH;col++){
        			final_image[rcvdRow][col] = partial_image[col];
        		}
  		
  		} while (count>0);
	}
	else{ //Slave code
		while(true){
			//Receive row from MASTER and check if it is the terminate tag
			MPI_Recv(&row, 1, MPI_INT,0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if(row==-1){
				break;
			}
			else{
				//Calculate greyscale values for this row and send it back to MASTER
				int i, j;
				c.imag = imag_min + (imag_max*row) /HEIGHT;
				for (i = 0; i < WIDTH; i++) {
				    
					c.real = real_min + (real_max*i) / WIDTH;
					partial_image[i] = cal_pixel(c);
				    
				}
				//Row tag is included to let MASTER know where to put the row of values in the final array
				MPI_Send(partial_image,WIDTH,MPI_INT,0,row,MPI_COMM_WORLD);
			}
		}
	}
	
	//End Timer
	double end_time = MPI_Wtime();
        
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank==0){
		double extime = (end_time - start_time);
		printf("the execution time is %f", extime);
		save_pgm("mandelbrot.pgm", final_image);
	}
	MPI_Finalize();

}	
