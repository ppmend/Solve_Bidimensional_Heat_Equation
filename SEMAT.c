//-----------------------------------------------------------------------------
/*
   Universidade Federal de Uberlândia
   Faculdade de Matemática
   Projeto de IC
Data: 18/08/2025
Autor: Pedro Paulo Nogueira Mendonça e Rafael Alves Figueiredo

Resolve a equação do calor 2D:

rho * Cp * psi_t = d( K * psi_x )/dx + d( K * psi_y )/dy + G

em um domínio [0,Lx] x [0, Ly]
usando o método de diferenças finitas e Gauss-Seidel
*/
//-----------------------------------------------------------------------------
// Bibliotecas do sistema

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

//-----------------------------------------------------------------------------

typedef struct {
  int **type;
  double **psi;
  double **K;
} StateField;

typedef struct {
  int Z;
  int Nx;
  int Ny;
  int Nt;
  double Lx;
  double Ly;
  double endtime;
  double T_infinity;
} params_t;

params_t *params; // Informação global

//-----------------------------------------------------------------------------
// Macros

#define PI acos(-1.0)

#define POINT_INVALID 0
#define POINT_VALID 1
#define POINT_DIRICHLET 2 // p = cte
#define POINT_NEUMAN1 4   // dp/dx = 0 || dp/dy = 0
#define POINT_NEUMAN2 8   // q_cond = q_conv

#define TOL 1.0e-8  // Tolerância do critério de parada da resolução do sistema linear
#define ITMAX 10.0e10   // Número máximo de iterações que o programa pode fazer
#define rho 8960.0   // Parâmetro da Equação do Calor
#define Cp 385.0     // Parâmetro da Equação do Calor
#define h 25.0      // Parâmetro da Equação do Calor

//-----------------------------------------------------------------------------

// Declaração de funções

StateField *SetDomain(int Z);
void Solver(StateField *field, int Z);
void write_vtk(StateField *field, int frame, int Z);
void printtype(int **type);
void printpsi(double **psi);
StateField *sf_alloc();
void freeField(StateField *field);
StateField *fieldCopy(StateField *field);
double CalcK(double T);
double calcG(double x);

//-----------------------------------------------------------------------------

int main() {

  // Declaração das variáveis
  //for(int Z = 26; Z <= 35; Z++ ){

  int Z = 50; // Parâmetro definidor do número de pontos do sistema

  params = calloc(1, sizeof(params_t)); // Nx, Ny, Lx, Ly são parâmetros globais

  StateField *field;

  field = SetDomain(Z);

  Solver(field, Z);

  printtype(field->type);
  printpsi(field->psi);

  freeField(field);
  free(params);

//}

  return 1;

}

//-----------------------------------------------------------------------------

StateField *SetDomain(int Z) { // Definindo domínio inicial

  StateField *field;
  params->Nx = 5 + 4 * Z ; // Nx = 5 + 4Z (Z é um inteiro)
  params->Ny = 6 + 5 * Z ; // Ny = 6 + 5Z (Z é um inteiro)
  params->Nt = 240;
  params->Lx = 50.0e-03; // Comprimento do sistema [m]
  params->Ly = 40.0e-03;
  params->endtime = 10.; // Tempo final [s]
  params->T_infinity = 298.15;

  field = sf_alloc();

  double x,y; // Coordenada dos pontos
  double dx,dy;

  dx = params->Lx / ( params->Nx-1.0 );
  dy = params->Ly / ( params->Ny-1.0 );

    for( int i = 0 ; i < params->Nx ; i++ ){ // Geometria que emula duas Aletas idênticas igualmente espaçadas e inicialmente a temperatura ambiente
    x = i*dx;

        for( int j = 0 ; j < params->Ny ; j++ ){
            y = j*dy;

            if( x < 0.25 * params->Lx && j != 0 && j != params->Ny - 1 && i != 0 ){

                field->type[i][j] = POINT_VALID;
                field->psi[i][j] = 298.15;
                field->K[i][j] = CalcK(field->psi[i][j]);

            }else if ( j == 0 && x < 0.25 * params->Lx ) {

                field->type[i][j] = POINT_NEUMAN1;
                field->psi[i][j] = 298.15;
                field->K[i][j] = CalcK(field->psi[i][j]);

            }else if( j == params->Ny - 1 && x < 0.25 * params->Lx ){

                field->type[i][j] = POINT_NEUMAN1;
                field->psi[i][j] = 298.15;
                field->K[i][j] = CalcK(field->psi[i][j]);

            }else if( i == 0 ){

                field->type[i][j] = POINT_NEUMAN1;
                field->psi[i][j] = 298.15;
                field->K[i][j] = CalcK(field->psi[i][j]);

            }else if( i == params->Nx - 1 && y >= 0.60 * params->Ly && y <= 0.80 * params->Ly){

                field->type[i][j] = POINT_NEUMAN2;
                field->psi[i][j] = 298.15;
                field->K[i][j] = CalcK(field->psi[i][j]);

            }else if( i == params->Nx - 1 && y >= 0.20 * params->Ly && y <= 0.40 * params->Ly){

                field->type[i][j] = POINT_NEUMAN2;
                field->psi[i][j] = 298.15;
                field->K[i][j] = CalcK(field->psi[i][j]);

            }else if( i == (params->Nx - 1)/4 && y >= 0. * params->Ly && y <= 0.20 * params->Ly){

                field->type[i][j] = POINT_NEUMAN2;
                field->psi[i][j] = 298.15;
                field->K[i][j] = CalcK(field->psi[i][j]);

            }else if( i == (params->Nx - 1)/4 && y >= 0.4 * params->Ly && y <= 0.60 * params->Ly){

                field->type[i][j] = POINT_NEUMAN2;
                field->psi[i][j] = 298.15;
                field->K[i][j] = CalcK(field->psi[i][j]);

            }else if( i == (params->Nx - 1)/4 && y >= 0.8 * params->Ly && y <= params->Ly){

                field->type[i][j] = POINT_NEUMAN2;
                field->psi[i][j] = 298.15;
                field->K[i][j] = CalcK(field->psi[i][j]);

            }else if( x >= 0.25 * params->Lx && j > (params->Ny - 1)/5  && j < 2 * (params->Ny - 1)/5 ){

                field->type[i][j] = POINT_VALID;
                field->psi[i][j] = 298.15;
                field->K[i][j] = CalcK(field->psi[i][j]);

            }else if( x >= 0.25 * params->Lx && j > 3 * (params->Ny - 1)/5  && j < 4 * (params->Ny - 1)/5){

                field->type[i][j] = POINT_VALID;
                field->psi[i][j] = 298.15;
                field->K[i][j] = CalcK(field->psi[i][j]);

            }else if( x >= 0.25 * params->Lx && j == (params->Ny - 1)/5 ){

                field->type[i][j] = POINT_NEUMAN2;
                field->psi[i][j] = 298.15;
                field->K[i][j] = CalcK(field->psi[i][j]);

            }else if( x >= 0.25 * params->Lx && j == 2 * (params->Ny - 1)/5){

                field->type[i][j] = POINT_NEUMAN2;
                field->psi[i][j] = 298.15;
                field->K[i][j] = CalcK(field->psi[i][j]);

            }else if( x >= 0.25 * params->Lx && j == 3 * (params->Ny - 1)/5){

                field->type[i][j] = POINT_NEUMAN2;
                field->psi[i][j] = 298.15;
                field->K[i][j] = CalcK(field->psi[i][j]);

            }else if( x >= 0.25 * params->Lx && j == 4 * (params->Ny - 1)/5){

                field->type[i][j] = POINT_NEUMAN2;
                field->psi[i][j] = 298.15;
                field->K[i][j] = CalcK(field->psi[i][j]);

            }else {

                field->type[i][j] = POINT_INVALID;
                field->psi[i][j] = 298.15;
                field->K[i][j] = CalcK(field->psi[i][j]); // No physical meaning
            }
        }
    }

  return field;
}

//-----------------------------------------------------------------------------

void Solver(StateField *field, int Z){

  StateField *field_old; // Campo temporário

  int k = 1;
  double x, y, t;
  double dx = (double) params->Lx / (params->Nx - 1.0);
  double dy = (double) params->Ly / (params->Ny - 1.0);
  double dt = (double) params->endtime / params->Nt;
  double part_sx = (double) dt / (rho * Cp * dx * dx);
  double part_sy = (double) dt / (rho * Cp * dy * dy);
  double psi_new, errorIt, sx, sy, wx, wy, rx, ry;
  int ItTotal = 0;

  write_vtk(field, 0, Z); // Salvar o campo inicial
  while(k <= params->Nt){

      field_old = fieldCopy(field); // Representa o tempo anterior (n <- n+1)
      t = (double) k * dt;

      for( int It = 0; It < ITMAX ; It++ ){
        errorIt = 0.;

        for( int i = 0 ; i < params->Nx ; i++ ){
          x = (double) i * dx;

          for( int j = 0 ; j < params->Ny ; j++ ){
            y = (double) j * dy;

            field->K[i][j] = CalcK(field->psi[i][j]);

            psi_new = field->psi[i][j];

            sx = field->K[i][j] * part_sx;
            sy = field->K[i][j] * part_sy;

            rx = ( 2 * dx * h ) / field->K[i][j];
            ry = ( 2 * dy * h ) / field->K[i][j];

            if( field->type[i][j] & POINT_VALID){

              wx = dt * (field->K[i+1][j] - field->K[i-1][j]) / (4. * dx * dx * rho * Cp);
              wy = dt * (field->K[i][j+1] - field->K[i][j-1]) / (4. * dy * dy * rho * Cp);

              psi_new = ( field_old->psi[i][j] + sx * (field->psi[i+1][j] + field->psi[i-1][j]) + sy * (field->psi[i][j+1] + field->psi[i][j-1]) + wx * ( field->psi[i+1][j] - field->psi[i-1][j] ) + wy * ( field->psi[i][j+1] - field->psi[i][j-1] ) + calcG(x) * dt / (rho * Cp) ) / ( 1. + 2. * sx + 2. * sy  );

              errorIt = fmax( errorIt , fabs(field->psi[i][j] - psi_new));
              field->psi[i][j] = psi_new;

            } else if( field->type[i][j] & POINT_NEUMAN1 ){ // Adiabatic conditions

              if( i == 0 ){ // Ghost point : [i-1][j]

                if( j == 0 ) {
                    psi_new = ( field_old->psi[i][j] + 2. * sx * field->psi[i+1][j] + 2 * sy * field->psi[i][j+1] + calcG(x) * dt / (rho * Cp) ) / ( 1. + 2. * sx + 2. * sy  );
                } else if (j == params->Ny - 1 ) {
                    psi_new = ( field_old->psi[i][j] + 2. * sx * field->psi[i+1][j] + 2 * sy * field->psi[i][j-1] + calcG(x) * dt / (rho * Cp) ) / ( 1. + 2. * sx + 2. * sy  );
                } else{
                    wy = dt * (field->K[i][j+1] - field->K[i][j-1]) / (4. * dy * dy * rho * Cp);
                    psi_new = ( field_old->psi[i][j] + 2. * sx * field->psi[i+1][j] + sy * (field->psi[i][j+1] + field->psi[i][j-1]) + wy * ( field->psi[i][j+1] - field->psi[i][j-1] ) + calcG(x) * dt / (rho * Cp) ) / ( 1. + 2. * sx + 2. * sy  );
                }

              } else if( j == 0 && i > 0){ // Ghost point : [i][j-1]

                wx = dt * (field->K[i+1][j] - field->K[i-1][j]) / (4. * dx * dx * rho * Cp);

                psi_new = ( field_old->psi[i][j] + sx * (field->psi[i+1][j] + field->psi[i-1][j]) + 2. * sy * field->psi[i][j+1] + wx * ( field->psi[i+1][j] - field->psi[i-1][j] ) + calcG(x) * dt / (rho * Cp) ) / ( 1. + 2. * sx + 2. * sy  );

              } else if( j == params->Ny - 1 && i > 0){ // Ghost point : [i][j+1]

                wx = dt * (field->K[i+1][j] - field->K[i-1][j]) / (4. * dx * dx * rho * Cp);

                psi_new = ( field_old->psi[i][j] + sx * (field->psi[i+1][j] + field->psi[i-1][j]) + 2. * sy * field->psi[i][j-1] + wx * ( field->psi[i+1][j] - field->psi[i-1][j] ) + calcG(x) * dt / (rho * Cp) ) / ( 1. + 2. * sx + 2. * sy  );

              }

            } else if( field->type[i][j] & POINT_NEUMAN2 ){ // q_cond = q_conv


              if(j == 0 && i == (params->Nx - 1) / 4 ){ // Ghost point : [i][j-1] and [i+1][j]

                wx = dt * (field->K[i][j] - field->K[i-1][j]) / (4. * dx * dx * rho * Cp);
                wy = dt * (field->K[i][j+1] - field->K[i][j]) / (4. * dy * dy * rho * Cp);

                psi_new = ( field_old->psi[i][j] + sy * (field->psi[i][j+1] + field->psi[i][j+1]) + sx * (2. * field->psi[i-1][j] - rx * params->T_infinity) - wx * rx * params->T_infinity + wy * ( field->psi[i][j+1] - field->psi[i][j+1] ) + calcG(x) * dt / (rho * Cp) ) / ( 1. + 2. * sx + 2. * sy - sx * rx - wx * rx );

              } else if(j == params->Ny - 1 && i == (params->Nx - 1) / 4 ){ // Ghost point : [i][j+1] and [i+1][j]

                wx = dt * (field->K[i][j] - field->K[i-1][j]) / (4. * dx * dx * rho * Cp);
                wy = dt * (field->K[i][j] - field->K[i][j-1]) / (4. * dy * dy * rho * Cp);

                psi_new = ( field_old->psi[i][j] + sy * (field->psi[i][j-1] + field->psi[i][j-1]) + sx * (2. * field->psi[i-1][j] - rx * params->T_infinity) - wx * rx * params->T_infinity + wy * ( field->psi[i][j-1] - field->psi[i][j-1] ) + calcG(x) * dt / (rho * Cp) ) / ( 1. + 2. * sx + 2. * sy - sx * rx - wx * rx );

              } else if(field->type[i][j+1] & POINT_INVALID){ // Ghost point : [i][j+1]

                wx = dt * (field->K[i+1][j] - field->K[i-1][j]) / (4. * dx * dx * rho * Cp);
                wy = dt * (field->K[i][j] - field->K[i][j-1]) / (4. * dy * dy * rho * Cp);

                psi_new = ( field_old->psi[i][j] + sx * (field->psi[i+1][j] + field->psi[i-1][j]) + sy * (2. * field->psi[i][j-1] - ry * params->T_infinity) - wy * ry * params->T_infinity + wx * ( field->psi[i+1][j] - field->psi[i-1][j] ) + calcG(x) * dt / (rho * Cp) ) / ( 1. + 2. * sx + 2. * sy - sy * ry - wy * ry );

              } else if(field->type[i][j-1] & POINT_INVALID){ // Ghost point : [i][j-1]

                wx = dt * (field->K[i+1][j] - field->K[i-1][j]) / (4. * dx * dx * rho * Cp);
                wy = dt * (field->K[i][j+1] - field->K[i][j]) / (4. * dy * dy * rho * Cp);

                psi_new = ( field_old->psi[i][j] + sx * (field->psi[i+1][j] + field->psi[i-1][j]) + sy * (2. * field->psi[i][j+1] + ry * params->T_infinity) + wy * ry * params->T_infinity + wx * ( field->psi[i+1][j] - field->psi[i-1][j] ) + calcG(x) * dt / (rho * Cp) ) / ( 1. + 2. * sx + 2. * sy + sy * ry + wy * ry );

              } else if( field->type[i+1][j] & POINT_INVALID || field->type[i-1][j] & POINT_VALID || field->type[i-1][j] & POINT_NEUMAN2 ){ // Ghost point : [i+1][j]

                wx = dt * (field->K[i][j] - field->K[i-1][j]) / (4. * dx * dx * rho * Cp);
                wy = dt * (field->K[i][j+1] - field->K[i][j-1]) / (4. * dy * dy * rho * Cp);

                psi_new = ( field_old->psi[i][j] + sy * (field->psi[i][j+1] + field->psi[i][j-1]) + sx * (2. * field->psi[i-1][j] - rx * params->T_infinity) - wx * rx * params->T_infinity + wy * ( field->psi[i][j+1] - field->psi[i][j-1] ) + calcG(x) * dt / (rho * Cp) ) / ( 1. + 2. * sx + 2. * sy - sx * rx - wx * rx );

              }
            }

            errorIt = fmax( errorIt , fabs(field->psi[i][j] - psi_new));
            field->psi[i][j] = psi_new;
          }
        }

        if( errorIt < TOL){
          printf("\nA convergencia no tempo t = %.2f (k = %d) foi alcancada apos %d iterações", t, k, It);
          ItTotal = ItTotal + It;
          break;
        }
      }

      freeField(field_old);

      write_vtk(field, k, Z);
      k = k + 1;
    }

  printf("\nSimulação z = %d finalizada com sucesso com %d iteracoes.\n", Z, ItTotal);

}

//-----------------------------------------------------------------------------

void write_vtk(StateField *field, int frame, int Z) {
  char filename[64];
  sprintf(filename, "z=%02doutput%04d.vtk", Z, frame);

  FILE *fp = fopen(filename, "w");
  if (!fp) {
    perror("Erro ao abrir arquivo");
    exit(1);
  }

  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "Heat equation solution\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET STRUCTURED_POINTS\n");
  fprintf(fp, "DIMENSIONS %d %d 1\n", params->Nx, params->Ny);
  fprintf(fp, "ORIGIN 0 0 0\n");
  fprintf(fp, "SPACING %lf %lf 1.0\n", params->Lx/(params->Nx-1), params->Ly/(params->Ny-1));
  fprintf(fp, "POINT_DATA %d\n", params->Nx * params->Ny);

  fprintf(fp, "SCALARS type float 1\n");
  fprintf(fp, "LOOKUP_TABLE default\n");
  for (int j = 0; j < params->Ny; j++) {
    for (int i = 0; i < params->Nx; i++) {
      fprintf(fp, "%i\n", field->type[i][j]);
    }
  }

  fprintf(fp, "SCALARS Temperature float 1\n");
  fprintf(fp, "LOOKUP_TABLE default\n");
  for (int j = 0; j < params->Ny; j++) {
    for (int i = 0; i < params->Nx; i++) {
      fprintf(fp, "%f\n", field->psi[i][j]);
    }
  }

  fclose(fp);
}

//-----------------------------------------------------------------------------

void printtype(int **type){

  printf("\n\n");
  for( int j = params->Ny - 1 ; j >= 0 ; j--){
	for( int i = 0 ; i < params->Nx ; i++ ){
	  printf("%d ", type[i][j]);
	}
    printf("\n");
  }

  return;
}

//-----------------------------------------------------------------------------

void printpsi(double **psi){

  printf("\n\n");
  for( int j = params->Ny - 1 ; j >= 0 ; j--){
	for( int i = 0 ; i < params->Nx ; i++ ){
	  printf("%.2f ", psi[i][j]); // Endtime print
	}
    printf("\n");
  }

  return;
}

//-----------------------------------------------------------------------------

StateField *sf_alloc() {

  StateField *field = calloc(1, sizeof(StateField));

  field->type = calloc(params->Nx, sizeof(int*));
  field->psi = calloc(params->Nx, sizeof(double*));
  field->K = calloc(params->Nx, sizeof(double*));

  for ( int i = 0 ; i < params->Nx ; i++ ) {
	field->type[i] = calloc(params->Ny, sizeof(int));
	field->psi[i] = calloc(params->Ny, sizeof(double));
	field->K[i] = calloc(params->Ny, sizeof(double));
  }

  return field;

}

//-----------------------------------------------------------------------------

void freeField(StateField *field) {

  for( int i = 0 ; i < params->Nx ; i++){
    free(field->psi[i]);
    free(field->type[i]);
    free(field->K[i]);
  }

  free(field->psi);
  free(field->type);
  free(field->K);

}

//-----------------------------------------------------------------------------

// Creates a copy of sf

StateField *fieldCopy(StateField *field) {

  StateField *field_new = calloc(1,sizeof(StateField));

  field_new = sf_alloc();

  for ( int i = 0 ; i < params->Nx ; i++ ) {
	memcpy(field_new->type[i], field->type[i], params->Ny * sizeof(int));
	memcpy(field_new->psi[i],  field->psi[i],  params->Ny * sizeof(double));
	memcpy(field_new->K[i],  field->K[i],  params->Ny * sizeof(double));
  }

return field_new;

}

//-----------------------------------------------------------------------------

double CalcK(double T){ // Condutividade [W/(m . K)]

    double a =  2227.664; // Bismuth PC Params and T = [Kelvin]
    double b = -0.627271;
    double c = 2.09554e-4;
    double d =  22.35452;

    //double a = 82.56648; // Cooper Params and T = [Kelvin]
    //double b = 0.262301;
    //double c = -4.06701e-4;
    //double d = 59.72934;

    return a * pow(T, b) * exp(c * T) * exp(d / T);
}

//-----------------------------------------------------------------------------

double calcG(double x){ // Termo fonte gerador de energia [W/m^3]

    double A = 1.0e9; // Param
    double B = 1.0e8; // Param
    double xf = 50.0e-3;

    double G = x * x * (B - A) / (xf * xf - 2. * xf) + 2 * x * (A - B) / (xf - 2.) + A;

    return G;

}

//-----------------------------------------------------------------------------
