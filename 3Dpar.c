/* Foudili Sirine
   Baaghbagha Hadjer
   Bio-Informatique 

   _____________________________________________
   Résolution de l'équation de chaleur
   _____________________________________________

   ---------------------
   Version parallèle
   ---------------------

   */






#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"


#define min(a,b) a <= b ? a : b


void computeNext(double*** x0, double*** x, int size_x, int size_y, int size_z, double dt, double hx, double hy, 
                 double hz, double* r, int me, int* xs, int* ys, int* zs, int* xe, int* ye, int* ze, double k0)
{

   int i, j, k;

   double diagx, diagy, diagz, weightx, weighty, weightz;

   double rk;



   diagx = -2.0 + hx*hx/(3*k0*dt);
   diagy = -2.0 + hy*hy/(3*k0*dt);
   diagz = -2.0 + hz*hz/(3*k0*dt);
   weightx = k0*dt/(hx*hx);
   weighty = k0*dt/(hy*hy);
   weightz = k0*dt/(hz*hz);

   for(k=zs[me];k<=ze[me];k++)
     for(j=ys[me];j<=ye[me];j++)
       for(i=xs[me];i<=xe[me];i++)
         x[i][j][k] = weightx *(x0[i-1][j][k] + x0[i+1][j][k] + x0[i][j][k]*diagx)
                    + weighty *(x0[i][j-1][k] + x0[i][j+1][k] + x0[i][j][k]*diagy)
                    + weightz *(x0[i][j][k-1] + x0[i][j][k+1] + x0[i][j][k]*diagz);


   *r = 0.0;
   for(k=zs[me];k<=ze[me];k++)
     for(j=ys[me];j<=ye[me];j++)
       for(i=xs[me];i<=xe[me];i++)
       {
        rk = x0[i][j][k] - x[i][j][k];
        *r  = *r + rk*rk;
        x0[i][j][k] = x[i][j][k];
       }
}    



void initValues(int nb_layers, double*** x0, int x_dim, int y_dim, int z_dim, double temp1_init, double temp2_init)
{

   int i, j, k, l;

   for(l=1;l<=nb_layers+1;l++)
   {
    for(j=(l-1);j<=y_dim-l;j++)
      for(k=(l-1);k<=z_dim-l;k++)
      {
       x0[l-1][j][k] = temp1_init;
       x0[x_dim-l][j][k] = temp1_init;
      }

    for(i=(l-1);i<=x_dim-l;i++)
      for(k=(l-1);k<=z_dim-l;k++)
      {
       x0[i][l-1][k] = temp1_init;
       x0[i][y_dim-l][k] = temp1_init;
      }

    for(i=(l-1);i<=x_dim-l;i++)
      for(j=(l-1);j<=y_dim-l;j++)
      {
       x0[i][j][l-1] = temp1_init;
       x0[i][j][z_dim-l] = temp1_init;
      }
   }

   for(i=nb_layers+1;i<=x_dim-(nb_layers+2);i++)
     for(j=nb_layers+1;j<=y_dim-(nb_layers+2);j++)
       for(k=nb_layers+1;k<=z_dim-(nb_layers+2);k++)
         x0[i][j][k] = temp2_init;
}


void processToMap(int me, int *xs, int *xe, int *ys, int *ye, int *zs, int *ze, int xcell, 
                  int ycell, int zcell, int x_domains, int y_domains, int z_domains)
{

   int i, j, k, l, m, p, v;

   for(i=0;i<=(z_domains*y_domains)-1;i++)
   {
    xs[i] = 2;
    xe[i] = xs[i]+xcell-1;
   }

   for(j=1;j<=x_domains-1;j++)
     for(k=0;k<=(z_domains*y_domains-1);k++)
     {
      xs[j*(z_domains*y_domains)+k] = xs[(j-1)*(z_domains*y_domains)]+xcell+2;
      xe[j*(z_domains*y_domains)+k] = xs[j*(z_domains*y_domains)]+xcell-1;
     }

   for(i=1;i<=y_domains;i++) {
      ys[(i-1)*z_domains] = y_domains*(ycell+2)-ycell*i-2*(i-1);
      ye[(i-1)*z_domains] = ys[(i-1)*z_domains]+ycell-1;

      for(l=1;l<=z_domains-1;l++) {
         ys[(i-1)*z_domains+l] = ys[(i-1)*z_domains];
         ye[(i-1)*z_domains+l] = ys[(i-1)*z_domains+l]+ycell-1;
      }
   }

   for(m=1;m<=y_domains;m++) {
      ys[(m-1)*z_domains] = y_domains*(ycell+2)-ycell*m-2*(m-1);
      ye[(m-1)*z_domains] = ys[(m-1)*z_domains]+ycell-1;

      for(i=1;i<=x_domains-1;i++) {
         ys[i*(y_domains*z_domains)+(m-1)*z_domains] = ys[(m-1)*z_domains];
         ye[i*(y_domains*z_domains)+(m-1)*z_domains] = ys[i*(y_domains*z_domains)+(m-1)*z_domains]+ycell-1;

         for(l=1;l<=z_domains-1;l++) {
            ys[i*(y_domains*z_domains)+(m-1)*z_domains+l] = ys[i*(y_domains*z_domains)+(m-1)*z_domains];
            ye[i*(y_domains*z_domains)+(m-1)*z_domains+l] = ys[i*(y_domains*z_domains)+(m-1)*z_domains+l]+ycell-1;
         }
      }
   }

   for(k=0;k<=y_domains-1;k++)
   {
    v = k*z_domains;
    zs[v] = 2;
    ze[v] = 2+zcell-1;
    for(p=1;p<=x_domains-1;p++)
    {
     zs[v+p*(y_domains*z_domains)] = zs[v];
     ze[v+p*(y_domains*z_domains)] = ze[v];
    }
   }

   for(m=1;m<=z_domains-1;m++)
      for(i=0;i<=y_domains-1;i++)
      {
       l = m+i*z_domains;
       zs[l] = zs[l-1]+zcell+2;
       ze[l] = zs[l]+zcell-1;
       for(v=1;v<=x_domains-1;v++)
       {
        zs[l+v*(y_domains*z_domains)] = zs[l];
        ze[l+v*(y_domains*z_domains)] = zs[l+v*(y_domains*z_domains)]+zcell-1;
       }
      }
}


int main(int argc, char *argv[])
{

   int size_x, size_y, size_z, me, x_domains, y_domains, z_domains;
   int size_x_glo, size_y_glo, size_z_glo;


   double ***x;
   double ***x0;
   double *x_all;
   double *x0_all;
   double *x_alloc;
   double *x0_alloc;
   double *xfinal;
   double *xtemp;

 
   int iconf[7];
   double conf[2];


   double dt, dt1, dt2, hx, hy, hz, min1, min2;


   double resLoc, result, epsilon;


   int convergence = 0;

   int i, j, k, l, p, v, m;


   double t;
   int step;

   int maxStep;

   double time_init, time_final;
   double elapsed_time;

   int nb_layers = 1;

 
   double temp1_init = 10.0;

   double temp2_init = -10.0;

   double k0 = 1;

   int sizes[3], subsizes1[3], subsizes2[3], subsizes3[3], starts[3];
   int nproc, ndims;
   MPI_Comm comm, comm3d;
   int dims[3];
   int periods[3];
   int reorganisation = 0;
   MPI_Datatype matrix_type_oxz, matrix_type_oxy, matrix_type_oyz;
   int S=0, E=1, N=2, W=3, Zd=4, Zu=5;
   int NeighBor[6];
   int xcell, ycell, zcell, size_tot_x, size_tot_y, size_tot_z;
   int *xs, *ys, *zs, *xe, *ye, *ze;

   MPI_Init(&argc, &argv);
   comm = MPI_COMM_WORLD;
   MPI_Comm_size(comm,&nproc);
   MPI_Comm_rank(comm,&me);

   if(me==0)
    readParam(iconf, conf);

   MPI_Bcast(iconf,7,MPI_INT,0,comm);
   MPI_Bcast(conf,2,MPI_DOUBLE,0,comm);

   size_x    = iconf[0];
   size_y    = iconf[1];
   size_z    = iconf[2];
   x_domains = iconf[3];
   y_domains = iconf[4];
   z_domains = iconf[5];
   maxStep   = iconf[6];
   dt1       = conf[0];
   epsilon   = conf[1];

   if((me==0) && (nproc!=(x_domains*y_domains*z_domains)))
    printf("Number of processes not equal to Number of subdomains\n");

   size_x_glo = size_x+2;
   size_y_glo = size_y+2;
   size_z_glo = size_z+2;
   hx = 1.0/(double)(size_x_glo);
   hy = 1.0/(double)(size_y_glo);
   hz = 1.0/(double)(size_z_glo);
   min1 = min(hx,hy);
   min2 = min(min1,hz);
   dt2  = 0.125*min2*min2*min2/k0;
   size_tot_x = size_x+2*x_domains+2;
   size_tot_y = size_y+2*y_domains+2;
   size_tot_z = size_z+2*z_domains+2;

   if(dt1>=dt2)
   {
    if(me==0)
    {
     printf("\n");
     printf("  Time step too large in 'param' file - Taking convergence criterion\n");
    }
    dt = dt2;
   }
   else dt = dt1;

   xfinal = malloc(size_x*size_y*size_z*sizeof(*xfinal));
   x_all = malloc(size_tot_x*size_tot_y*size_tot_z*sizeof(*x_all)); 
   x0_all = malloc(size_tot_x*size_tot_y*size_tot_z*sizeof(*x0_all));
   x_alloc = x_all;
   x0_alloc = x0_all;

   x = malloc(size_tot_x*sizeof(*x));
   x0 = malloc(size_tot_x*sizeof(*x0));

   for (i=0;i<size_tot_x;i++)
   {
    x[i] = malloc(size_tot_y*sizeof(**x));
    x0[i] = malloc(size_tot_y*sizeof(**x0));
    for(j=0;j<size_tot_y;j++)
    {
     x[i][j] = x_alloc;
     x0[i][j] = x0_alloc;
     x_alloc += size_tot_z;
     x0_alloc += size_tot_z;
    }
   }

   
   xs = malloc(nproc*sizeof(int));
   xe = malloc(nproc*sizeof(int));
   ys = malloc(nproc*sizeof(int));
   ye = malloc(nproc*sizeof(int));
   zs = malloc(nproc*sizeof(int));
   ze = malloc(nproc*sizeof(int));

   periods[0] = 0;
   periods[1] = 0;
   periods[2] = 0;

   ndims = 3;
   dims[0] = x_domains;
   dims[1] = y_domains;
   dims[2] = z_domains;

   MPI_Cart_create(comm, ndims, dims, periods, reorganisation, &comm3d);

   NeighBor[0] = MPI_PROC_NULL;
   NeighBor[1] = MPI_PROC_NULL;
   NeighBor[2] = MPI_PROC_NULL;
   NeighBor[3] = MPI_PROC_NULL;
   NeighBor[4] = MPI_PROC_NULL;
   NeighBor[5] = MPI_PROC_NULL;

   MPI_Cart_shift(comm3d, 0, 1, &NeighBor[W], &NeighBor[E]);

   MPI_Cart_shift(comm3d, 1, 1, &NeighBor[S], &NeighBor[N]);

   MPI_Cart_shift(comm3d, 2, 1, &NeighBor[Zd], &NeighBor[Zu]);
   xcell = (size_x/x_domains);
   ycell = (size_y/y_domains);
   zcell = (size_z/z_domains);

   xtemp = malloc(xcell*ycell*zcell*sizeof(*xtemp));

   processToMap(me, xs, xe, ys, ye, zs, ze, xcell, ycell, zcell, x_domains, y_domains, z_domains);

   sizes[0] = size_tot_x;
   sizes[1] = size_tot_y;
   sizes[2] = size_tot_z;

   starts[0] = 0;
   starts[1] = 0;
   starts[2] = 0;

   subsizes1[0] = xcell;
   subsizes1[1] = 1;
   subsizes1[2] = zcell;

   MPI_Type_create_subarray(3, sizes, subsizes1, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix_type_oxz);
   MPI_Type_commit(&matrix_type_oxz);

   subsizes2[0] = 1;
   subsizes2[1] = ycell;
   subsizes2[2] = zcell;

   MPI_Type_create_subarray(3, sizes, subsizes2, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix_type_oyz);
   MPI_Type_commit(&matrix_type_oyz);

   subsizes3[0] = xcell;
   subsizes3[1] = ycell;
   subsizes3[2] = 1;

   MPI_Type_create_subarray(3, sizes, subsizes3, starts, MPI_ORDER_C, MPI_DOUBLE, &matrix_type_oxy);
   MPI_Type_commit(&matrix_type_oxy);

   initValues(nb_layers, x0, size_tot_x, size_tot_y, size_tot_z, temp1_init, temp2_init);

   updateBound(x0, size_tot_x, size_tot_y, size_tot_z, NeighBor, comm3d,
               matrix_type_oxz, matrix_type_oxy, matrix_type_oyz, me, xs, ys, zs, xe, ye, ze);

   step = 0;
   t = 0.0;

   time_init = MPI_Wtime();

   while(!convergence)
   {  
      step = step + 1;
      t = t + dt ;
      computeNext(x0, x, size_tot_x, size_tot_y, size_tot_z, dt, hx, hy, hz, &resLoc, me, xs, ys, zs, xe, ye, ze, k0);

      updateBound(x0, size_tot_x, size_tot_y, size_tot_z, NeighBor, comm3d, 
                  matrix_type_oxz, matrix_type_oxy, matrix_type_oyz, me, xs, ys, zs, xe, ye, ze);

      MPI_Allreduce(&resLoc, &result, 1, MPI_DOUBLE, MPI_SUM, comm);

      result = sqrt(result);

      if ((result<epsilon) || (step>maxStep)) break;
   }

   i = 1;
   for(k=zs[me];k<=ze[me];k++)
   {
    l = 1;
    for(j=ys[me];j<=ye[me];j++)
    {
     for(m=0;m<=xcell-1;m++)
       xtemp[(l-1)*xcell+(i-1)*xcell*ycell+m] = x0[xs[me]+m][j][k];
     l = l+1;
    }
    i = i+1;
   }
   MPI_Gather(xtemp, xcell*ycell*zcell, MPI_DOUBLE, xfinal, xcell*ycell*zcell, MPI_DOUBLE, 0, comm);

   time_final = MPI_Wtime();
   elapsed_time = time_final - time_init;

   if(me == 0)
    printf("  Temps d'execution = %15.6f\n",elapsed_time);


   /* Free arrays */
   for (i=0;i<=size_tot_x-1;i++)
   {
    free(x[i]);
    free(x0[i]);
   }
   
   free(x);
   free(x0);
   free(x_all);
   free(x0_all);
   free(xfinal);
   free(xtemp);
   free(xs);
   free(xe);
   free(ys);
   free(ye);
   free(zs);
   free(ze);

   /* Free matrices type */
   MPI_Type_free(&matrix_type_oxz);
   MPI_Type_free(&matrix_type_oxy);
   MPI_Type_free(&matrix_type_oyz);

   MPI_Finalize();

   return 0;
}