
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"
#include <math.h>

// Define output file name
#define OUTPUT_FILE "stencil.pgm"
#define MASTER 0

void stencil(unsigned short nx, unsigned short ny, float * restrict  image, float * restrict  tmp_image);
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

  int ranknx = ceil(ny/nnodes);
  int rankny = nx;
  int nxWhole = ranknx + 2;

  if(myrank == 0){
    nxWhole -= 1;
  }

  if(myrank == nnodes -1 ){
    nxWhole = nx - ((nnodes-1)*rankny)+1;
  }

  // Allocate the image
  float * restrict image = malloc(sizeof(float)*nx*ny);
  float * restrict tmp_image = malloc(sizeof(float)*nx*ny);

  float * restrict rankImg = malloc(sizeof(float)*rankny*nxWhole);
  float * restrict tmp_rankImg = malloc(sizeof(float)*rankny*nxWhole);

  init_image(nx, ny, image, tmp_image);

  if(myrank == 0){
    for(int i = 0; i < rankny*nxWhole; i++){
      rankImg[i] = image[i];
      tmp_rankImg[i] = image[i];
    }
  }else{
    for(int i = 0; i < rankny*nxWhole; i++){
      rankImg[i] = image[(myrank*rankny*ranknx-rankny)+i];
      tmp_rankImg[i] = image[(myrank*rankny*ranknx-rankny)+i];
    }
  }

  if(myrank == nnodes-1){
    ranknx = nx - ((nnodes-1)*ranknx);
  }

////////////////////////////////////////
  double tic = wtime();
  for (unsigned short t = 0; t < niters; ++t) {

    stencil(nxWhole, rankny, rankImg, tmp_rankImg);

    if(myrank > 0){
      MPI_Sendrecv(&tmp_rankImg[rankny], rankny, MPI_FLOAT, myrank-1, tag,
                    &tmp_rankImg[0], rankny, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, &status);
    }

    if(myrank < nnodes-1){
      MPI_Sendrecv(&tmp_rankImg[rankny*nxWhole-(2*rankny)], rankny, MPI_FLOAT, myrank+1, tag,
                    &tmp_rankImg[rankny*nxWhole], rankny, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, &status);
    }
/////////////////////////////////////////////////
    stencil(nxWhole, rankny, tmp_rankImg, rankImg);

    if(myrank > 0){
      MPI_Sendrecv(&rankImg[rankny], rankny, MPI_FLOAT, myrank-1, tag,
                    &rankImg[0], rankny, MPI_FLOAT, myrank-1, tag, MPI_COMM_WORLD, &status);
    }

    if(myrank < nnodes-1){
      MPI_Sendrecv(&rankImg[rankny*nxWhole-(2*rankny)], rankny, MPI_FLOAT, myrank+1, tag,
                    &rankImg[rankny*nxWhole], rankny, MPI_FLOAT, myrank+1, tag, MPI_COMM_WORLD, &status);
    }
  }
  double toc = wtime();

  if(myrank == 0){
    for(int i = 0; i < rankny*ranknx; i++){
      image[i] = rankImg[i];
    }

    for(int n = 1; n < nnodes; n++){
      for(int x = 0; x < ranknx; x++){
        MPI_Recv(&image[n*(rankny*ranknx)+ x*rankny], rankny, MPI_FLOAT, n, tag, MPI_COMM_WORLD, &status);
      }
    }
  }else{
    for(int x = 1; x < nxWhole; x++){
      MPI_Send(&rankImg[x*rankny], rankny, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
    }
  }
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

  void stencil(unsigned short nx, unsigned short ny, float * restrict  image, float * restrict  tmp_image) {

    //if(rank == 0){
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
    //}

    //if(rank == nnodes-1){
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
    //}



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
