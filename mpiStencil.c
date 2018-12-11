
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"
#include <math.h>

// Define output file name
#define OUTPUT_FILE "stencil.pgm"
#define MASTER 0

//void stencil(unsigned short nx, unsigned short ny, float * restrict image, float * restrict tmp_image);
void stencil(unsigned short nx, unsigned short ny, float * restrict  image, float * restrict  tmp_image, int rank, int nnodes);
void init_image(unsigned short nx, unsigned short ny, float * restrict  image, float * restrict  tmp_image);
void output_image(const char * file_name, unsigned short nx, unsigned short ny, float * restrict  image);
double wtime(void);

int main(int argc, char *argv[]) {

  //setup MPI
  int myrank;
  int nnodes;
  int tag = 0;
  MPI_Status status;

  MPI_Init( &argc, &argv );

  MPI_Comm_size( MPI_COMM_WORLD, &nnodes );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  // Check usage
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  // Initiliase problem dimensions from command line arguments
  int nx = atoi(argv[1]);
  float ny = atoi(argv[2]);
  int niters = atoi(argv[3]);

  int ranknx = nx;
  int rankny = ceil(ny/nnodes);
  int nyWhole = rankny + 2;
  //tesing MPI
  int rankStart, rankEnd;

  if(myrank == 0){
    nyWhole -= 1;
    rankStart = 0;
  }else{
    rankStart = myrank * (rankny*ranknx) - ranknx;
  }

  if(myrank == nnodes-1){
    nyWhole -= 1;
    rankEnd = nx*ny-1;
    //rankny = ((rankEnd+1) - rankStart) / ranknx;
  }else{
    rankEnd = ((myrank+1) * (ranknx*rankny)) - 1 + ranknx;
  }

  // Allocate the image
  float * restrict image = malloc(sizeof(float)*nx*ny);
  float * restrict tmp_image = malloc(sizeof(float)*nx*ny);

  float * restrict rankImg = malloc(sizeof(float)*ranknx*nyWhole);
  float * restrict tmp_rankImg = malloc(sizeof(float)*ranknx*nyWhole);

  printf("Rank = %d, Rank Start = %d, Rank End = %d, nx = %d, ny = %d, Allocated mem = %d\n"
          , myrank, rankStart, rankEnd, ranknx, nyWhole, ranknx*nyWhole);

  // Set the input image
  init_image(nx, ny, image, tmp_image);

  //set up each ranks section of the image
  if(myrank == 0){
    for(int i = 0; i < ranknx*nyWhole; i++){
      rankImg[i] = image[i];
      tmp_rankImg[i] = image[i];
    }
  }else{
    for(int i = 0; i < ranknx*nyWhole; i++){
      rankImg[i] = image[myrank*ranknx*rankny-ranknx+i];
      tmp_rankImg[i] = image[myrank*ranknx*rankny-ranknx+i];
    }
  }

  int even = myrank%2 == 0 ? 1 : 0;   //check if rank is even

  // Call the stencil kernel
  double tic = wtime();
  for (unsigned short t = 0; t < niters; ++t) {
    //stencil(nx, ny, image, tmp_image);
    //stencil(nx, ny, tmp_image, image);
    stencil(ranknx, nyWhole, rankImg, tmp_rankImg, myrank, nnodes);

    //if(even){
    if(myrank != 0){
      MPI_Sendrecv(&tmp_rankImg[ranknx], ranknx, MPI_FLOAT, myrank-1, tag,
                    &tmp_rankImg[0], ranknx, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, &status);
    }


      // if(myrank != 0){
      //   MPI_Send(&tmp_rankImg[ranknx], ranknx, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD);
      //   /*for(int i = ranknx; i < 2*ranknx; i++){
      //     // printf("1) Rank %d sent %0.1lf to rank %d\n", myrank, rankImg[i], myrank-1 );
      //     MPI_Send(&rankImg[i], 1, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD);
      //   }*/
      // }
      // if(myrank != 0){
      //   MPI_Recv(&tmp_rankImg[0], ranknx, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, &status);
      //   /*for(int i = 0; i < ranknx; i++){
      //     // printf("4) Rank %d recieved %0.1lf from rank %d\n", myrank, rankImg[i], myrank-1);
      //     MPI_Recv(&rankImg[i], 1, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, &status);
      //   }*/
      // }
      if(myrank != nnodes-1){
        MPI_Sendrecv(&tmp_rankImg[ranknx*nyWhole-(2*ranknx)], ranknx, MPI_FLOAT, myrank+1, tag,
                      &tmp_rankImg[ranknx*nyWhole-ranknx], ranknx, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, &status);
      }


      // if(myrank != nnodes-1){
      //   MPI_Send(&tmp_rankImg[ranknx*nyWhole-(2*ranknx)], ranknx, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD);
      //   /*for(int i = ranknx*nyWhole-2*ranknx; i < ranknx*nyWhole-ranknx; i++){
      //     // printf("2) Rank %d sent %0.1lf to rank %d\n", myrank, rankImg[i], myrank+1);
      //     MPI_Send(&rankImg[i], 1, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD);
      //   }*/
      // }
      // if(myrank != nnodes-1){
      //   MPI_Recv(&tmp_rankImg[ranknx*nyWhole-ranknx], ranknx, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, &status);
      //   /*for(int i = ranknx*nyWhole-ranknx; i < ranknx*nyWhole; i++){
      //     // printf("3) Rank %d recieved %0.1lf from rank %d\n", myrank, rankImg[i], myrank+1);
      //     MPI_Recv(&rankImg[i], 1, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, &status);
      //   }*/
      // }


    //}else{
    if(myrank != nnodes-1){
      MPI_Sendrecv(&tmp_rankImg[ranknx*nyWhole-(2*ranknx)], ranknx, MPI_FLOAT, myrank+1, tag,
                    &tmp_rankImg[ranknx*nyWhole-ranknx], ranknx, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, &status);
    }


      // if(myrank != nnodes-1){
      //   MPI_Recv(&tmp_rankImg[ranknx*nyWhole-ranknx], ranknx, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, &status);
      //   /*for(int i = ranknx*nyWhole-ranknx; i < ranknx*nyWhole; i++){
      //     // printf("5) Rank %d recieved %0.1lf from rank %d\n", myrank, rankImg[i], myrank+1);
      //     MPI_Recv(&rankImg[i], 1, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, &status);
      //   }*/
      // }
      // if(myrank != nnodes-1){
      //   MPI_Send(&tmp_rankImg[ranknx*nyWhole-(2*ranknx)], ranknx, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD);
      //   /*for(int i = ranknx*nyWhole-2*ranknx; i < ranknx*nyWhole-ranknx; i++){
      //     //printf("8) Rank %d sent %0.1lf to rank %d\n", myrank, rankImg[i], myrank+1);
      //     MPI_Send(&rankImg[i], 1, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD);
      //   }*/
      // }
      if(myrank != 0){
        MPI_Sendrecv(&tmp_rankImg[ranknx], ranknx, MPI_FLOAT, myrank-1, tag,
                      &tmp_rankImg[0], ranknx, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, &status);
      }


      // MPI_Recv(&tmp_rankImg[0], ranknx, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, &status);
      // /*for(int i = 0; i < ranknx; i++){
      //   // printf("6) Rank %d recieved %0.1lf from rank %d\n", myrank, rankImg[i], myrank-1);
      //   MPI_Recv(&rankImg[i], 1, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, &status);
      // }*/
      //
      // MPI_Send(&tmp_rankImg[ranknx], ranknx, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD);
      // /*for(int i = ranknx; i < 2*ranknx; i++){
      //   //printf("7) Rank %d sent %0.1lf to rank %d\n", myrank, rankImg[i], myrank-1 );
      //   MPI_Send(&rankImg[i], 1, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD);
      // }*/
    //}
    stencil(ranknx, nyWhole, tmp_rankImg, rankImg, myrank, nnodes);
////////////////////////////////////////////////////////////////////////////////////////////

//if(even){
if(myrank != 0){
  MPI_Sendrecv(&tmp_rankImg[ranknx], ranknx, MPI_FLOAT, myrank-1, tag,
                &tmp_rankImg[0], ranknx, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, &status);
}


  // if(myrank != 0){
  //   MPI_Send(&tmp_rankImg[ranknx], ranknx, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD);
  //   /*for(int i = ranknx; i < 2*ranknx; i++){
  //     // printf("1) Rank %d sent %0.1lf to rank %d\n", myrank, rankImg[i], myrank-1 );
  //     MPI_Send(&rankImg[i], 1, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD);
  //   }*/
  // }
  // if(myrank != 0){
  //   MPI_Recv(&tmp_rankImg[0], ranknx, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, &status);
  //   /*for(int i = 0; i < ranknx; i++){
  //     // printf("4) Rank %d recieved %0.1lf from rank %d\n", myrank, rankImg[i], myrank-1);
  //     MPI_Recv(&rankImg[i], 1, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, &status);
  //   }*/
  // }
  if(myrank != nnodes-1){
    MPI_Sendrecv(&tmp_rankImg[ranknx*nyWhole-(2*ranknx)], ranknx, MPI_FLOAT, myrank+1, tag,
                  &tmp_rankImg[ranknx*nyWhole-ranknx], ranknx, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, &status);
  }


  // if(myrank != nnodes-1){
  //   MPI_Send(&tmp_rankImg[ranknx*nyWhole-(2*ranknx)], ranknx, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD);
  //   /*for(int i = ranknx*nyWhole-2*ranknx; i < ranknx*nyWhole-ranknx; i++){
  //     // printf("2) Rank %d sent %0.1lf to rank %d\n", myrank, rankImg[i], myrank+1);
  //     MPI_Send(&rankImg[i], 1, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD);
  //   }*/
  // }
  // if(myrank != nnodes-1){
  //   MPI_Recv(&tmp_rankImg[ranknx*nyWhole-ranknx], ranknx, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, &status);
  //   /*for(int i = ranknx*nyWhole-ranknx; i < ranknx*nyWhole; i++){
  //     // printf("3) Rank %d recieved %0.1lf from rank %d\n", myrank, rankImg[i], myrank+1);
  //     MPI_Recv(&rankImg[i], 1, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, &status);
  //   }*/
  // }


//}else{
if(myrank != nnodes-1){
  MPI_Sendrecv(&tmp_rankImg[ranknx*nyWhole-(2*ranknx)], ranknx, MPI_FLOAT, myrank+1, tag,
                &tmp_rankImg[ranknx*nyWhole-ranknx], ranknx, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, &status);
}


  // if(myrank != nnodes-1){
  //   MPI_Recv(&tmp_rankImg[ranknx*nyWhole-ranknx], ranknx, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, &status);
  //   /*for(int i = ranknx*nyWhole-ranknx; i < ranknx*nyWhole; i++){
  //     // printf("5) Rank %d recieved %0.1lf from rank %d\n", myrank, rankImg[i], myrank+1);
  //     MPI_Recv(&rankImg[i], 1, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, &status);
  //   }*/
  // }
  // if(myrank != nnodes-1){
  //   MPI_Send(&tmp_rankImg[ranknx*nyWhole-(2*ranknx)], ranknx, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD);
  //   /*for(int i = ranknx*nyWhole-2*ranknx; i < ranknx*nyWhole-ranknx; i++){
  //     //printf("8) Rank %d sent %0.1lf to rank %d\n", myrank, rankImg[i], myrank+1);
  //     MPI_Send(&rankImg[i], 1, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD);
  //   }*/
  // }
  if(myrank != 0){
    MPI_Sendrecv(&tmp_rankImg[ranknx], ranknx, MPI_FLOAT, myrank-1, tag,
                  &tmp_rankImg[0], ranknx, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, &status);
  }


  // MPI_Recv(&tmp_rankImg[0], ranknx, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, &status);
  // /*for(int i = 0; i < ranknx; i++){
  //   // printf("6) Rank %d recieved %0.1lf from rank %d\n", myrank, rankImg[i], myrank-1);
  //   MPI_Recv(&rankImg[i], 1, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, &status);
  // }*/
  //
  // MPI_Send(&tmp_rankImg[ranknx], ranknx, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD);
  // /*for(int i = ranknx; i < 2*ranknx; i++){
  //   //printf("7) Rank %d sent %0.1lf to rank %d\n", myrank, rankImg[i], myrank-1 );
  //   MPI_Send(&rankImg[i], 1, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD);
  // }*/
//}
  }

  if(myrank == 0){
    for(int i = 0; i < ranknx*rankny; i++){
      image[i] = rankImg[i];
    }

    for(int n = 1; n < nnodes; n++){
      for(int y = 0; y < rankny; y++){
        MPI_Recv(&image[n*(ranknx*rankny) + y*ranknx], ranknx, MPI_FLOAT, n, tag, MPI_COMM_WORLD, &status);
      }
    }
  }else{
    for(int y = 1; y < nyWhole; y++){
      MPI_Send(&rankImg[y*ranknx], ranknx, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
    }
  }
  double toc = wtime();

  printf("Rank %d image stiched, time = %lf\n", myrank, toc-tic);

  if(myrank == 0){
    // Output
    printf("------------------------------------\n");
    printf(" runtime: %lf s\n", toc-tic);
    printf("------------------------------------\n");

    output_image(OUTPUT_FILE, nx, ny, image);
  }

  free(image);
  free(rankImg);
  free(tmp_rankImg);

  MPI_Finalize();
}

void stencil(unsigned short nx, unsigned short ny, float * restrict  image, float * restrict  tmp_image, int rank, int nnodes) {

    //topLeft
    tmp_image[0] = image[0] * 0.6f;                     //order of each of these sections could matter with respect to caching
    tmp_image[0] += image[1]* 0.1f;
    tmp_image[0] += image[nx]* 0.1f;

    //topRight
    tmp_image[nx-1] = image[nx-1] * 0.6f;
    tmp_image[nx-1] += image[nx-2]* 0.1f;
    tmp_image[nx-1] += image[(2*nx)-1]* 0.1f;

    for(int topEdge = 1; topEdge < nx-1; ++topEdge){
      tmp_image[topEdge] = image[topEdge] * 0.6f;
      tmp_image[topEdge] += image[topEdge + 1]* 0.1f;
      tmp_image[topEdge] += image[topEdge - 1]* 0.1f;
      tmp_image[topEdge] += image[topEdge + nx]*0.1f;
    }

    //bottomLeft
    tmp_image[nx*(ny-1)] = image[nx*(ny-1)] * 0.6f;
    tmp_image[nx*(ny-1)] += image[nx*(ny-1)+1]* 0.1f;
    tmp_image[nx*(ny-1)] += image[nx*(ny-2)]* 0.1f;

    //bottomRight
    tmp_image[(nx*ny)-1] = image[(nx*ny)-1] * 0.6f;
    tmp_image[(nx*ny)-1] += image[(nx*ny)-2]* 0.1f;
    tmp_image[(nx*ny)-1] += image[(nx*ny)-(nx+1)]* 0.1f;

    for(int bottomEdge = 1; bottomEdge < nx-1; ++bottomEdge){
      tmp_image[nx*(ny-1)+bottomEdge] = image[nx*(ny-1)+bottomEdge] * 0.6f;
      tmp_image[nx*(ny-1)+bottomEdge] += image[(nx*(ny-1)+bottomEdge) + 1]* 0.1f;
      tmp_image[nx*(ny-1)+bottomEdge] += image[(nx*(ny-1)+bottomEdge) - 1]* 0.1f;
      tmp_image[nx*(ny-1)+bottomEdge] += image[(nx*(ny-1)+bottomEdge) - nx]*0.1f;
    }

  for(int leftEdge = 1; leftEdge < ny-1; ++leftEdge){
    tmp_image[nx*leftEdge] = image[nx*leftEdge] * 0.6f;
    tmp_image[nx*leftEdge] += image[nx*leftEdge + 1]* 0.1f;
    tmp_image[nx*leftEdge] += image[nx*leftEdge + nx]* 0.1f;
    tmp_image[nx*leftEdge] += image[nx*leftEdge - nx]*0.1f;
  }

  //unsigned short rightPixel = nx+(nx-1);
  for(int rightEdge = 1; rightEdge < ny-1; ++rightEdge){
    tmp_image[rightEdge*nx+(nx-1)] = image[rightEdge*nx+(nx-1)] * 0.6f;
    tmp_image[rightEdge*nx+(nx-1)] += image[(rightEdge*nx+(nx-1)) - 1]* 0.1f;
    tmp_image[rightEdge*nx+(nx-1)] += image[(rightEdge*nx+(nx-1)) + nx]* 0.1f;
    tmp_image[rightEdge*nx+(nx-1)] += image[(rightEdge*nx+(nx-1)) - nx]*0.1f;
  }

  for (int i = 1; i < ny-1; ++i) {
    for (int j = 1; j < nx-1; ++j) {

      tmp_image[j+i*nx] = image[j+i*nx] * 0.6f;
      tmp_image[j+i*nx] += image[(j+i*nx)+1] * 0.1f;
      tmp_image[j+i*nx] += image[(j+i*nx) - 1] * 0.1f;
      tmp_image[j+i*nx] += image[(j+i*nx) + nx] * 0.1f;
      tmp_image[j+i*nx] += image[(j+i*nx) - nx] * 0.1f;
    }
  }
}

// Create the input image
void init_image(unsigned short nx, unsigned short ny, float * restrict  image, float * restrict  tmp_image) {
  // Zero everything
  for (unsigned short j = 0; j < ny; ++j) {
    for (unsigned short i = 0; i < nx; ++i) {
      image[j+i*ny] = 0.0f;
      tmp_image[j+i*ny] = 0.0f;
    }
  }

  // Checkerboard
  for (unsigned short j = 0; j < 8; ++j) {
    for (unsigned short i = 0; i < 8; ++i) {
      for (unsigned short jj = j*ny/8; jj < (j+1)*ny/8; ++jj) {
        for (unsigned short ii = i*nx/8; ii < (i+1)*nx/8; ++ii) {
          if ((i+j)%2)
          image[jj+ii*ny] = 100.0f;
        }
      }
    }
  }
}

// Routine to output the image in Netpbm grayscale binary image format
void output_image(const char * file_name, unsigned short nx, unsigned short ny, float * restrict  image) {

  // Open output file
  FILE *fp = fopen(file_name, "w");
  if (!fp) {
    fprintf(stderr, "Error: Could not open %s\n", OUTPUT_FILE);
    exit(EXIT_FAILURE);
  }

  // Ouptut image header
  fprintf(fp, "P5 %d %d 255\n", nx, ny);

  // Calculate maximum value of image
  // This is used to rescale the values
  // to a range of 0-255 for output
  float maximum = 0.0f;
  for (unsigned short j = 0; j < ny; ++j) {
    for (unsigned short i = 0; i < nx; ++i) {
      if (image[j+i*ny] > maximum)
        maximum = image[j+i*ny];
    }
  }

  // Output image, converting to numbers 0-255
  for (unsigned short j = 0; j < ny; ++j) {
    for (unsigned short i = 0; i < nx; ++i) {
      fputc((char)(255.0f*image[j+i*ny]/maximum), fp);
    }
  }

  // Close the file
  fclose(fp);

}

// Get the current time in seconds since the Epoch
double wtime(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec*1e-6;
}
