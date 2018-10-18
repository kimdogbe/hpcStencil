
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

// Define output file name
#define OUTPUT_FILE "stencil.pgm"

void stencil(unsigned short nx, unsigned short ny, float * restrict  image, float * restrict  tmp_image);
void init_image(unsigned short nx, unsigned short ny, float * restrict  image, float * restrict  tmp_image);
void output_image(const char * file_name, unsigned short nx, unsigned short ny, float * restrict  image);
double wtime(void);

int main(int argc, char *argv[]) {

  // Check usage
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  // Initiliase problem dimensions from command line arguments
  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  int niters = atoi(argv[3]);

  // Allocate the image
  float *image = malloc(sizeof(float)*nx*ny);
  float *tmp_image = malloc(sizeof(float)*nx*ny);

  // Set the input image
  init_image(nx, ny, image, tmp_image);

  // Call the stencil kernel
  double tic = wtime();
  for (int t = 0; t < niters; ++t) {
    stencil(nx, ny, image, tmp_image);
    stencil(nx, ny, tmp_image, image);
  }
  double toc = wtime();


  // Output
  printf("------------------------------------\n");
  printf(" runtime: %lf s\n", toc-tic);
  printf("------------------------------------\n");

  output_image(OUTPUT_FILE, nx, ny, image);
  free(image);
}

void stencil(unsigned short nx, unsigned short ny, float * restrict  image, float * restrict  tmp_image) {

  //topLeft
  tmp_image[0] = image[0] * 0.6f;                     //order of each of these sections could matter with respect to caching
  tmp_image[0] += image[1]* 0.1f;
  tmp_image[0] += image[nx]* 0.1f;

  //topRight
  tmp_image[nx-1] = image[nx-1] * 0.6f;
  tmp_image[nx-1] += image[nx-2]* 0.1f;
  tmp_image[nx-1] += image[(nx-1)+nx]* 0.1f;

  //bottomLeft
  tmp_image[nx*(nx-1)] = image[nx*(nx-1)] * 0.6f;
  tmp_image[nx*(nx-1)] += image[nx*(nx-1)+1]* 0.1f;
  tmp_image[nx*(nx-1)] += image[nx*(nx-1)-nx]* 0.1f;

  //bottomRight
  tmp_image[(nx*nx)-1] = image[(nx*nx)-1] * 0.6f;
  tmp_image[(nx*nx)-1] += image[(nx*nx)-2]* 0.1f;
  tmp_image[(nx*nx)-1] += image[(nx*nx)-nx]* 0.1f;

  for(unsigned short topEdge = 1; topEdge < nx-1; topEdge++){
    tmp_image[topEdge] = image[topEdge] * 0.6f;
    tmp_image[topEdge] += image[topEdge + nx]* 0.1f;
    tmp_image[topEdge] += image[topEdge - 1]* 0.1f;
    tmp_image[topEdge] += image[topEdge + 1]*0.1f;
  }

  for(unsigned short bottomEdge = 1; bottomEdge < nx-1; bottomEdge++){
    tmp_image[nx*(nx-1)+bottomEdge] = image[nx*(nx-1)+bottomEdge] * 0.6f;
    tmp_image[nx*(nx-1)+bottomEdge] += image[nx*(nx-1)+bottomEdge + 1]* 0.1f;
    tmp_image[nx*(nx-1)+bottomEdge] += image[nx*(nx-1)+bottomEdge - 1]* 0.1f;
    tmp_image[nx*(nx-1)+bottomEdge] += image[nx*(nx-1)+bottomEdge - nx]*0.1f;
  }

  for(unsigned short leftEdge = 1; leftEdge < nx-1; leftEdge++){
    tmp_image[nx*leftEdge] = image[nx*leftEdge] * 0.6f;
    tmp_image[nx*leftEdge] += image[nx*leftEdge + 1]* 0.1f;
    tmp_image[nx*leftEdge] += image[nx*leftEdge - nx]* 0.1f;
    tmp_image[nx*leftEdge] += image[nx*leftEdge + nx]*0.1f;
  }

  //unsigned short rightPixel = nx+(nx-1);
  for(unsigned short rightEdge = 1; rightEdge < nx-1; rightEdge++){
    tmp_image[rightEdge*nx+(nx-1)] = image[rightEdge*nx+(nx-1)] * 0.6f;
    tmp_image[rightEdge*nx+(nx-1)] += image[rightEdge*nx+(nx-1) - 1]* 0.1f;
    tmp_image[rightEdge*nx+(nx-1)] += image[rightEdge*nx+(nx-1) + nx]* 0.1f;
    tmp_image[rightEdge*nx+(nx-1)] += image[rightEdge*nx+(nx-1) - nx]*0.1f;
  }

  for (unsigned short i = 1; i < nx-1; ++i) {
    for (unsigned short j = 1; j < nx-1; ++j) {

      tmp_image[j+i*nx] = image[j+i*nx] * 0.6f;
      tmp_image[j+i*nx] += image[j  +(i-1)*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j  +(i+1)*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j-1+i*nx] * 0.1f;
      tmp_image[j+i*nx] += image[j+1+i*nx] * 0.1f;
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
  for (int j = 0; j < 8; ++j) {
    for (int i = 0; i < 8; ++i) {
      for (int jj = j*ny/8; jj < (j+1)*ny/8; ++jj) {
        for (int ii = i*nx/8; ii < (i+1)*nx/8; ++ii) {
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
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      if (image[j+i*ny] > maximum)
        maximum = image[j+i*ny];
    }
  }

  // Output image, converting to numbers 0-255
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
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
