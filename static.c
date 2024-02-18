#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>

#define WIDTH 640
#define HEIGHT 480
#define MAX_ITER 255

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
	
	//Calculate size of partition that each processor will take
	int partition = HEIGHT/(size);
	
	//Two arrays, a local one for the different partitions and final_image for the master
	
	int final_image[HEIGHT][WIDTH];
	int partial_image[partition][WIDTH];
	struct complex c;
	
	//Start Timer
	//Iterate through the rows of respective processor
	//Fill partial_image array with the mandelbrot set values (0-255)
	double start_time = MPI_Wtime();
        
        
	//Initialize boundaries for each processor
	//Process i will start from the row (i*partition)
	int start = rank*partition;
	int end = start + partition;
	int i, j;
	for (i = 0; i < partition; i++) {
	    for (j = 0; j < WIDTH; j++) {
		c.real = -1.5 + (3.0*j) / WIDTH; //Applying the appropriate mapping from 
		c.imag = -1.7 + (3.0*(i + start)) /HEIGHT; //the display values to coordinates on the complex plane
		partial_image[i][j] = cal_pixel(c);
	    }
	}
	
	MPI_Gather(partial_image,partition*WIDTH,MPI_INT,final_image, partition*WIDTH,MPI_INT,0,MPI_COMM_WORLD);
        
        //End Timer	
	double end_time = MPI_Wtime();
	
	

	if(rank==0){
		double extime = (end_time - start_time);
		printf("the execution time is %f from rank %d\n", extime,rank);
		save_pgm("mandelbrot.pgm", final_image);
	}
	
	MPI_Finalize();

}	
