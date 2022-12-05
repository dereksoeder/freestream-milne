#include <stddef.h>
#include <unistd.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include "FileIO.h"
#include "FSConfig.h"
#include <math.h>

void writeScalarToFile(float *var, char name[255], const parameters & params)
{
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  float DX = params.DX;
  float DY = params.DY;
  float DETA = params.DETA;
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (size_t ix = 0; ix < DIM_X; ix++)
  {
    for (size_t iy = 0; iy < DIM_Y; iy++)
    {
      for (size_t ieta = 0; ieta < DIM_ETA; ieta++)
      {
        float x = (float)ix * DX  - (((float)(DIM_X-1)) / 2.0 * DX);
        x = DX * roundf(x / DX);
        float y = (float)iy * DY  - (((float)(DIM_Y-1)) / 2.0 * DY);
        y = DY * roundf(y / DY);
        float eta = (float)ieta * DETA  - (((float)(DIM_ETA-1)) / 2.0 * DETA);
        eta = DETA * roundf(eta / DETA);

        size_t is = ix + (DIM_X * iy) + (DIM_X * DIM_Y * ieta); //the column packed index spanning x, y, z

        myfile << x << " " << y << " " << eta << " " << var[is] << "\n";
      }
    }
  }
  myfile.close();
}

void writeVectorToFile(float **var, char name[255], size_t idx, const parameters & params)
{
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  float DX = params.DX;
  float DY = params.DY;
  float DETA = params.DETA;
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (size_t ix = 0; ix < DIM_X; ix++)
  {
    for (size_t iy = 0; iy < DIM_Y; iy++)
    {
      for (size_t ieta = 0; ieta < DIM_ETA; ieta++)
      {
        float x = (float)ix * DX  - (((float)(DIM_X-1)) / 2.0 * DX);
        x = DX * roundf(x / DX); //rounding for regularly spaced values
        float y = (float)iy * DY  - (((float)(DIM_Y-1)) / 2.0 * DY);
        y = DY * roundf(y / DY);
        float eta = (float)ieta * DETA  - (((float)(DIM_ETA-1)) / 2.0 * DETA);
        eta = DETA * roundf(eta / DETA);

        size_t is = ix + (DIM_X * iy) + (DIM_X * DIM_Y * ieta); //the column packed index spanning x, y, z

        myfile << x << " " << y << " " << eta << " " << var[idx][is] << "\n";
      }
    }
  }
  myfile.close();
}

//this function writes the transverse density of a variable at z = 0
// as regularly spaced values
void writeScalarToFileProjection(float *var, char name[255], const parameters & params)
{
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (size_t iy = 0; iy < DIM_Y; iy++)
  {
    for (size_t ix = 0; ix < DIM_X; ix++)
    {
      size_t ieta = (DIM_ETA - 1) / 2; // at eta = 0
      size_t is = ix + (DIM_X * iy) + (DIM_X * DIM_Y * ieta); //the column packed index spanning x, y, eta
      myfile << var[is] << " "; //different columns for x values
    }
    myfile << "\n"; // different rows correspond to different y values
  }
  myfile.close();
}

void writeVectorToFileProjection(float **var, char name[255], size_t idx, const parameters & params)
{
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (size_t iy = 0; iy < DIM_Y; iy++)
  {
    for (size_t ix = 0; ix < DIM_X; ix++)
    {
      size_t ieta = (DIM_ETA - 1) / 2; //at eta = 0
      size_t is = ix + (DIM_X * iy) + (DIM_X * DIM_Y * ieta); //the column packed index spanning x, y, z
      myfile << var[idx][is] << " "; //different columns for x values
    }
    myfile << "\n"; // different rows correspond to different y values
  }
  myfile.close();
}

void readDensityFile(float *density, char name[255], const parameters & params)
{
  size_t DIM_X = params.DIM_X;
  size_t DIM_Y = params.DIM_Y;
  size_t DIM_ETA = params.DIM_ETA;
  float DX = params.DX;
  float DY = params.DY;
  float DETA = params.DETA;
  float xmin = (-1.0) * ((float)(DIM_X-1) / 2.0) * DX;
  float ymin = (-1.0) * ((float)(DIM_Y-1) / 2.0) * DY;
  float etamin = (-1.0) * ((float)(DIM_ETA-1) / 2.0) * DETA;
  float x, y, eta, value;

  char filename[255] = "";
  sprintf(filename, "%s.dat", name);
  std::ifstream infile;
  infile.open(filename);
  if (!infile)
  {
    printf("Couldn't open initial profile!\n");
    exit(1);
  }
  while (infile >> x >> y >> eta >> value)
  {
    size_t ix = (size_t)round((x - xmin) / DX);
    size_t iy = (size_t)round((y - ymin) / DY);
    size_t ieta = (size_t)round((eta - etamin) / DETA);
    if ((ix >= DIM_X) || (iy >= DIM_Y) || (ieta >= DIM_ETA))  // basic check for over-/underflow
    {
      printf("Coordinates (%f, %f, %f) out of bounds (%zu, %zu, %zu)!\n", x, y, eta, params.DIM_X, params.DIM_Y, params.DIM_ETA);
      exit(1);
    }
    size_t is = ix + (DIM_X * iy) + (DIM_X * DIM_Y * ieta);
    density[is] = value;
  }
  infile.close();
}

void readInParameters(struct parameters &params, const char * inputPath)
{
  char dummyChar[255];
  int dummyInt;
  //long dummyLong;
  float dummyFloat;
  size_t dummySize = 0;

  if (inputPath == nullptr)
  {
    inputPath = "freestream_input";
  }

  FILE *fileIn;
  fileIn = fopen(inputPath,"r");

  if (fileIn == nullptr)
  {
    printf("Couldn't open '%s' to read parameters. Using default values!\n", inputPath);
  }

  else
  {
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.OUTPUTFORMAT = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.BARYON = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.IC_ENERGY = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.IC_BARYON = dummyInt;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.ETA_WIDTH = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.ETA_FLAT = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.SIGMA = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.SIGMA_B = dummyFloat;
    fscanf(fileIn, "%s\t%zu\n", dummyChar, &dummySize);
    params.DIM_X = dummySize;
    fscanf(fileIn, "%s\t%zu\n", dummyChar, &dummySize);
    params.DIM_Y = dummySize;
    fscanf(fileIn, "%s\t%zu\n", dummyChar, &dummySize);
    params.DIM_ETA = dummySize;
    fscanf(fileIn, "%s\t%zu\n", dummyChar, &dummySize);
    params.DIM_RAP = dummySize;
    fscanf(fileIn, "%s\t%zu\n", dummyChar, &dummySize);
    params.DIM_PHIP = dummySize;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.DX = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.DY = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.DETA = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.DRAP = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.DTAU = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.TAU0 = dummyFloat;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.EOS_TYPE = dummyInt;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.E_FREEZE = dummyFloat;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.VISCOUS_MATCHING = dummyInt;

    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.E_DEP_FS = dummyInt;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.E_R = dummyFloat;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.TAU_R = dummyFloat;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.ALPHA = dummyFloat;

    fclose(fileIn);
  }

}
/*
void readEoSTable()
{

}
*/
