/* Various utility functions for dealing with basis sets */


/* Copyright (c) 2007 MJ Rutter 
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License version 2
 * as published by the Free Software Foundation.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA  02110-1301, USA.
 */ 


#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "cellsym.h"


/* Horrible hack to pass from main() here the fact that the user
 * cannot tell his left from his right, and wants the world to
 * remain that way. */
extern int lhs_fudge;

/* Check whether set is right- or left-handed */
int is_rhs(double b[3][3]){
  if ((b[0][0]*(b[1][1]*b[2][2]-b[1][2]*b[2][1])
      -b[0][1]*(b[1][0]*b[2][2]-b[1][2]*b[2][0])
      +b[0][2]*(b[1][0]*b[2][1]-b[1][1]*b[2][0]))<0) return(0);
  else return(1);
}

/* Update global reciprocal basis set to reflect global real basis set */
void real2rec(void){
  int i,j,k;
  double v;

  v=0.0;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      v+=basis[i][j]*basis[i][j];

  if (v==0.0) error_exit("all cell axes are length zero: no basis found?");

  for(i=0;i<3;i++){
    j=(i+1)%3;
    k=(i+2)%3;
    recip[i][0]=basis[j][1]*basis[k][2]-basis[j][2]*basis[k][1];
    recip[i][1]=-basis[j][0]*basis[k][2]+basis[j][2]*basis[k][0];
    recip[i][2]=basis[j][0]*basis[k][1]-basis[j][1]*basis[k][0];
  }

  v=0.0;
  for(i=0;i<3;i++) v+=basis[0][i]*recip[0][i];

  if (v==0.0) error_exit("unit cell volume is zero.");

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      recip[i][j]=recip[i][j]/v;

  cell_vol=v;
}

/* Update global fractional atomic co-ordinates to reflect global
 * absolute atomic co-ordinates. Assumes global reciprocal basis
 * is correct
 */
void addfrac(void){
  int i,j,k;

  for(i=0;i<natoms;i++){
    for(j=0;j<3;j++){
      atoms[i].frac[j]=0;
      for(k=0;k<3;k++)
         atoms[i].frac[j]+=atoms[i].abs[k]*recip[j][k];
    }
  }
}

/* Update global absolute atomic co-ordinates to reflect global
 * fractional atomic co-ordinates
 */
void addabs(void){
  int i,j,k;

  for(i=0;i<natoms;i++){
    for(j=0;j<3;j++){
      atoms[i].abs[j]=0;
      for(k=0;k<3;k++)
         atoms[i].abs[j]+=atoms[i].frac[k]*basis[k][j];
    }
  }
}

/* Fix global co-ordinates to force everything into 1st unit cell
 * Assume all global cell descriptions are consistent
 */
void reduce_cell(void){
  int i,j,fixed=0;

  for(i=0;i<natoms;i++){
    for(j=0;j<3;j++){
      if ((atoms[i].frac[j]>=0)&&(atoms[i].frac[j]<=1)) continue;
      fixed=1;
      atoms[i].frac[j]=fmod(atoms[i].frac[j],1.0);
      if (atoms[i].frac[j]<0) atoms[i].frac[j]+=1.0;
    }
  }

  if (fixed) addabs();
}

void abc2cart(double *abc, double basis[3][3]){
/* convert from a,b,c,alpha,beta,gamma to Cartesian basis */
  double alpha,beta,gamma,x;
  int i;

  if (debug>2)fprintf(stderr,"abc2cart: original basis:\n%f %f %f\n%f %f %f\n",
     abc[0],abc[1],abc[2],abc[3],abc[4],abc[5]);

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

  alpha=abc[3]*M_PI/180;
  beta=abc[4]*M_PI/180;
  gamma=abc[5]*M_PI/180;

/* a lies along x axis */

  basis[0][0]=abc[0];
  basis[0][1]=0.0;
  basis[0][2]=0.0;

/* b is in xy plane and angle gamma to x */

  basis[1][0]=abc[1]*cos(gamma);
  basis[1][1]=abc[1]*sin(gamma);
  basis[1][2]=0.0;

/* a,b,c is a right-hand set */

  basis[2][0]=abc[2]*cos(beta);

/* basis[2].basis[1] = abc[2] abc[1] cos(alpha) */

  x=abc[2]*abc[1]*cos(alpha)-basis[1][0]*basis[2][0];
  basis[2][1]=0;
  if (fabs(x)>1e-20){
    if (fabs(basis[1][1])>1e-30){
      basis[2][1]=x/basis[1][1];
    }else{
      error_exit("impossible problem in abc2cart, perhaps gamma is zero.\n");
    }
  }

/* And mod(basis[2][])=abc[2] */

  basis[2][2]=sqrt(abc[2]*abc[2]-basis[2][0]*basis[2][0]
                                -basis[2][1]*basis[2][1]);

/* Worry about handedness */

  real2rec();
  if (cell_vol<0){
    for(i=0;i<3;i++) basis[2][i]*=-1;
    cell_vol*=-1;
  }

  if(debug>2){
    int i;
    fprintf(stderr,"abc2cart: final basis:\n");
      for(i=0;i<=2;i++)
        fprintf(stderr,"%f %f %f\n",basis[i][0],basis[i][1],basis[i][2]);
  }

}

/* Convert from cartesian basis set to a,b,c,alpha,beta,gamma */
void cart2abc(double basis[3][3], double *abc, int fix){
  int i,j,k;
  double dtmp;

  if (debug>2) fprintf(stderr,"cart2abc called with fix=%d\n",fix);

  for(i=0;i<3;i++)
    abc[i]=sqrt(basis[i][0]*basis[i][0]+basis[i][1]*basis[i][1]+
                basis[i][2]*basis[i][2]);

  for(i=3;i<6;i++){
    j=(i+1)%3;
    k=(i+2)%3;
    abc[i]=acos((basis[j][0]*basis[k][0]+basis[j][1]*basis[k][1]+
                 basis[j][2]*basis[k][2])/(abc[j]*abc[k]))*180/M_PI;
  }

  /* We may now have a different orientation to the original, so: */

  if (fix){
    if (is_rhs(basis)==0){ /* We are trying to convert a lhs of vectors
                            * to abc format.
               * We should warn people,
               * Exchange second and third axes
               * Exchange second and third fractional co-ords
               */
      if (lhs_fudge) {
         fprintf(stderr,"Need rh set of vectors for abc output.\n"
                        "Have lh set, and exchange prohibitted...\n");
      }
      else {
        if (debug) fprintf(stderr,"Need rh set of vectors for abc output. "
                           "Exchanging basis vectors %d and %d.\n",
                           (1+preserve_c)%3+1,(2+preserve_c)%3+1);

        real2rec();

        for(i=0;i<3;i++){
          dtmp=basis[(1+preserve_c)%3][i];
          basis[(1+preserve_c)%3][i]=basis[(2+preserve_c)%3][i];
          basis[(2+preserve_c)%3][i]=dtmp;
        }
        real2rec();
        for(i=0;i<natoms;i++){
          dtmp=atoms[i].frac[(1+preserve_c)%3];
          atoms[i].frac[(1+preserve_c)%3]=atoms[i].frac[(2+preserve_c)%3];
          atoms[i].frac[(2+preserve_c)%3]=dtmp;
        }

        dtmp=abc[(1+preserve_c)%3];
        abc[(1+preserve_c)%3]=abc[(2+preserve_c)%3];
        abc[(2+preserve_c)%3]=dtmp;

        dtmp=abc[(1+preserve_c)%3+3];
        abc[(1+preserve_c)%3+3]=abc[(2+preserve_c)%3+3];
        abc[(2+preserve_c)%3+3]=dtmp;

      }
    } /* if (fix) */

    abc2cart(abc,basis);
    real2rec();
    addabs();
  }

}


/* Convert from cartesian basis set to a,b,c,alpha,beta,gamma */
void basis2abc(double basis[3][3], double abc[6]){
  int i,j,k;

  for(i=0;i<3;i++)
    abc[i]=sqrt(basis[i][0]*basis[i][0]+basis[i][1]*basis[i][1]+
                basis[i][2]*basis[i][2]);

  for(i=3;i<6;i++){
    j=(i+1)%3;
    k=(i+2)%3;
    abc[i]=acos((basis[j][0]*basis[k][0]+basis[j][1]*basis[k][1]+
                 basis[j][2]*basis[k][2])/(abc[j]*abc[k]))*180/M_PI;
  }
}
