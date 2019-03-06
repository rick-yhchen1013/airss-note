
/* Copyright (c) 2013 MJ Rutter 
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

/* routine to identify a symmetry operation */

void vcross(double a[3],double b[3],double c[3]);
int minvert(double m[3][3]);
double vmod2(double v[3]);

void ident_sym(struct sym_op *s,FILE *out){
  int i,j,inv,mult,screw;
  double m[3][3];
  double a[3],af[3],off[3];
  double v[3],v2[3];
  double det,angle,x;

  inv=0;
  screw=0;
  for(i=0;i<3;i++) off[i]=0;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      m[i][j]=s->mat[i][j];

  /* Calculate determinant */

  vcross(m[1],m[2],v);
  det=m[0][0]*v[0]+m[0][1]*v[1]+m[0][2]*v[2];

  if(!aeq(fabs(det),1)){
    fprintf(stderr,"Surprise: determinant is %f\n",det);
    return;
  }

  if (aeq(det,1)){
    if (debug>2) fprintf(stderr,"Looks like a rotation matrix\n");
  }

  if (aeq(det,-1)){
    if (debug>2) fprintf(stderr,"Looks like an inverse rotation matrix\n");
    inv=1;
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	m[i][j]*=-1;
  }

  angle=acos(0.5*(m[0][0]+m[1][1]+m[2][2]-1));
  angle*=180/3.1415926;

  if (debug>2) fprintf(stderr,"Angle looks like %f\n",angle);

  mult=0;
  if(aeq(angle,0)) mult=1;
  if(aeq(angle,180)) mult=2;
  if(aeq(angle,120)) mult=3;
  if(aeq(angle,90)) mult=4;
  if(aeq(angle,60)) mult=6;
  if (mult==0){
    fprintf(stderr,"Impossible angle %f in ident_sym\n",angle);
    return;
  }

  /* Now need axis. It will be e-vector with e-value +1 */

  if (mult==1){
    v[0]=v[1]=0;
    v[2]=1;
  }
  else{
    for(i=0;i<3;i++)
      m[i][i]-=1;
    vcross(m[0],m[1],v);
    if (aeq(vmod2(v),0)) vcross(m[1],m[2],v);
    if (aeq(vmod2(v),0)) vcross(m[2],m[0],v);
    if (aeq(vmod2(v),0)){
      fprintf(stderr,"Surprise! 1\n");
      return;
    }
  }

  for(i=0;i<3;i++) a[i]=v[i];
  if (debug>2) fprintf(stderr,"Axis looks like (%f,%f,%f)\n",a[0],a[1],a[2]);

  /* Try that in a fractional basis */

  for(i=0;i<3;i++)
    v2[i]=v[0]*recip[i][0]+v[1]*recip[i][1]+v[2]*recip[i][2];
  x=1e10;
  if ((fabs(v2[0])>tol)&&(fabs(v2[0])<x)) x=fabs(v2[0]);
  if ((fabs(v2[1])>tol)&&(fabs(v2[1])<x)) x=fabs(v2[1]);
  if ((fabs(v2[2])>tol)&&(fabs(v2[2])<x)) x=fabs(v2[2]);
  for(i=0;i<3;i++) v2[i]/=x;
  for(i=0;i<3;i++) a[i]/=x;
    
  for(i=0;i<3;i++) af[i]=v2[i];
  if (debug>2)
    fprintf(stderr,"Axis looks like (%f,%f,%f)\n",af[0],af[1],af[2]);

  /* Now deal with any offset */

  if (s->tr){
    for(i=0;i<3;i++) v[i]=s->tr[i];
    if(!(aeq(v[0],0)&&aeq(v[1],0)&&aeq(v[2],0))){
      x=v[0]*a[0]+v[1]*a[1]+v[2]*a[2];
      if((inv==0)&&(!aeq(x,0))){   /* we have a screw */
	x/=vmod2(a);  /* dot with unit vector parallel to axis, then divide by
		       * repeat length in that direction. Hence no sqrt. */
	for(i=0;i<3;i++) v[i]-=x*a[i];
	x=fmod(x,1);
	if (x<0) x+=1;
	if (x>1-tol) x=0;
	x*=mult;
	if (!aeq(x,floor(x+tol)))
	  fprintf(stderr,"Screw problem, x=%f\n",x);
	else
	  screw=x+tol;
      }
    }
    if(!(aeq(v[0],0)&&aeq(v[1],0)&&aeq(v[2],0))){ /* we have a translation */
      //      fprintf(stderr,"Have tr=(%f,%f,%f)\n",v[0],v[1],v[2]);
      if((inv==1)&&(aeq(angle,0))){ /* An inversion point */
	for(i=0;i<3;i++) off[i]=0.5*v[i];
      }
      if((inv==1)&&(aeq(angle,180))){ /* A mirror -- two evals are one */
	for(i=0;i<3;i++) off[i]=0.5*v[i];  /* Hmmm. Should take only bit
					    * parallel to a */
      }
      else{
	for(i=0;i<3;i++)
	  for(j=0;j<3;j++)
	    m[i][j]=-s->mat[i][j];
	for(i=0;i<3;i++) m[i][i]+=1.0;
	/* We now have 1-R, and the translation is we have is the result
         * of (1-R)t, where t is the translation we want. Unfortunately,
         * if R is a proper rotation, it has an evalue of 1, so 1-R
         * has a zero evalue and is uninvertable. Fix this.by adding a matrix
	 * with two zero evals, and one non-zero eval whose evec is along
	 * the rotation axis. This squashes the zero evalue, and leaves
	 * the others and their evectors alone.
	 */
	if (inv==0){
	  for(i=0;i<3;i++)
	    for(j=0;j<3;j++)
	      m[i][j]+=a[i]*a[j];
#if 0
	  if(minvert(m)){
	    for(j=0;j<3;j++) m[0][j]-=a[j];
	    if(minvert(m)){
	      for(j=0;j<3;j++) m[1][j]-=a[j];
	      if(minvert(m)){
		fprintf(stderr,"Failed to construct invertable matrix.\n");
		fprintf(stderr,"Seem to have a %d axis along (%f,%f,%f)"
			" inv=%d\n",mult,
			af[0],af[1],af[2],inv);
		return;
	      }
	    }
	  }
#endif
	}

	if(minvert(m)){
	  fprintf(stderr,"Surprise! Can't invert matrix\n");
	  fprintf(stderr,"Seem to have a %d axis along (%f,%f,%f)"
		  " inv=%d\n",mult,
		  af[0],af[1],af[2],inv);
	  return;
	}
	for(i=0;i<3;i++)
	  off[i]=m[i][0]*v[0]+m[i][1]*v[1]+m[i][2]*v[2];
      }
    }
  }

  if (inv) mult*=-1;
  for(i=0;i<3;i++)
    v[i]=off[0]*recip[i][0]+off[1]*recip[i][1]+off[2]*recip[i][2];

  if (screw)
    fprintf(out,"%d_%d axis along (%6.3f,%6.3f,%6.3f)",mult,screw,
	    af[0],af[1],af[2]);
  else
    fprintf(out,"%2d  axis along (%6.3f,%6.3f,%6.3f)",mult,
	    af[0],af[1],af[2]);

  if (aeq(v[0],0)&&aeq(v[1],0)&&aeq(v[2],0))
    fprintf(out," through (0,0,0)\n");
  else
    fprintf(out," through (%f,%f,%f)\n",v[0],v[1],v[2]);
}

void mpr(double m[3][3]){
  int i;
  for(i=0;i<3;i++)
    fprintf(stderr,"[%7f %7f %7f]\n",m[i][0],m[i][1],m[i][2]);
}

void vcross(double a[3],double b[3],double c[3]){
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-b[2]*a[0];
  c[2]=a[0]*b[1]-b[0]*a[1];
}

double vmod2(double v[3]){
  return (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

int minvert(double m[3][3]){
  int i,j;
  double det;
  double v[3],c[3][3];
#if 0
  double t[3][3];
#endif

  vcross(m[1],m[2],v);
  det=m[0][0]*v[0]+m[0][1]*v[1]+m[0][2]*v[2];
  if (aeq(det,0)) {
#if 0
    if (debug){
      fprintf(stderr,"Singular matrix:\n");
      mpr(m);
    }
#endif
    return 1;
  }
  /* Transpose of cofactor matrix */
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c[j][i]=m[(1+i)%3][(1+j)%3]*m[(2+i)%3][(2+j)%3]-
	m[(2+i)%3][(1+j)%3]*m[(1+i)%3][(2+j)%3];

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c[i][j]=c[i][j]/det;

#if 0
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      t[i][j]=c[i][0]*m[0][j]+c[i][1]*m[1][j]+c[i][2]*m[2][j];

  fprintf(stderr,"Identify in minv:\n");
  mpr(t);
#endif

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      m[i][j]=c[i][j];

  return 0;
}
