
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
  float * restrict image = malloc(sizeof(float)*nx*ny);
  float * restrict tmp_image = malloc(sizeof(float)*nx*ny);

  // Set the input image
  init_image(nx, ny, image, tmp_image);

  // Call the stencil kernel
  double tic = wtime();
  for (unsigned short t = 0; t < niters; ++t) {
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
  tmp_image[0] += image[ny]* 0.1f;

  //topRight
  tmp_image[ny-1] = image[ny-1] * 0.6f;
  tmp_image[ny-1] += image[ny-2]* 0.1f;
  tmp_image[ny-1] += image[(2*ny)-1]* 0.1f;

  //bottomLeft
  tmp_image[ny*(nx-1)] = image[ny*(nx-1)] * 0.6f;
  tmp_image[ny*(nx-1)] += image[ny*(nx-1)+1]* 0.1f;
  tmp_image[ny*(nx-1)] += image[ny*(nx-2)]* 0.1f;

  //bottomRight
  tmp_image[(ny*nx)-1] = image[(ny*nx)-1] * 0.6f;
  tmp_image[(ny*nx)-1] += image[(ny*nx)-2]* 0.1f;
  tmp_image[(ny*nx)-1] += image[(ny*nx)-(ny+1)]* 0.1f;

  for(int leftEdge = 1; leftEdge < nx-1; ++leftEdge){
    tmp_image[ny*leftEdge] = image[ny*leftEdge] * 0.6f;
    tmp_image[ny*leftEdge] += image[ny*leftEdge + 1]* 0.1f;
    tmp_image[ny*leftEdge] += image[ny*leftEdge + ny]* 0.1f;
    tmp_image[ny*leftEdge] += image[ny*leftEdge - ny]*0.1f;
  }

  //unsigned short rightPixel = nx+(nx-1);
  for(int rightEdge = 1; rightEdge < nx-1; ++rightEdge){
    tmp_image[rightEdge*ny+(ny-1)] = image[rightEdge*ny+(ny-1)] * 0.6f;
    tmp_image[rightEdge*ny+(ny-1)] += image[(rightEdge*ny+(ny-1)) - 1]* 0.1f;
    tmp_image[rightEdge*ny+(ny-1)] += image[(rightEdge*ny+(ny-1)) + nx]* 0.1f;
    tmp_image[rightEdge*ny+(ny-1)] += image[(rightEdge*ny+(ny-1)) - nx]*0.1f;
  }

  for(int topEdge = 1; topEdge < ny-1; ++topEdge){
    tmp_image[topEdge] = image[topEdge] * 0.6f;
    tmp_image[topEdge] += image[topEdge + 1]* 0.1f;
    tmp_image[topEdge] += image[topEdge - 1]* 0.1f;
    tmp_image[topEdge] += image[topEdge + ny]*0.1f;
  }

  for(int bottomEdge = 1; bottomEdge < ny-1; ++bottomEdge){
    tmp_image[ny*(nx-1)+bottomEdge] = image[ny*(nx-1)+bottomEdge] * 0.6f;
    tmp_image[ny*(nx-1)+bottomEdge] += image[(ny*(nx-1)+bottomEdge) + 1]* 0.1f;
    tmp_image[ny*(nx-1)+bottomEdge] += image[(ny*(nx-1)+bottomEdge) - 1]* 0.1f;
    tmp_image[ny*(nx-1)+bottomEdge] += image[(ny*(nx-1)+bottomEdge) - ny]*0.1f;
  }


  for (int i = 1; i < nx-1; ++i) {
    for (int j = 1; j < ny-1; ++j) {

      tmp_image[j+i*ny] = image[j+i*ny] * 0.6f;
      tmp_image[j+i*ny] += image[(j+i*ny)+1] * 0.1f;
      tmp_image[j+i*ny] += image[(j+i*ny) - 1] * 0.1f;
      tmp_image[j+i*ny] += image[(j+i*ny) + ny] * 0.1f;
      tmp_image[j+i*ny] += image[(j+i*ny) - ny] * 0.1f;
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
