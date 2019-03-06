/* Write a PDB file */


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

#include "cellsym.h"

void ident_sym(struct sym_op *s, FILE *out);

void cell_write_kpts(FILE* outfile);

void cell_write(FILE* outfile){
  int i;

  if (title)
    fprintf(outfile,"#TITL %s\n",title);

  fprintf(outfile,"%%block LATTICE_CART\nang\n");

  for(i=0;i<3;i++)
    fprintf(outfile,"%f %f %f\n",
                  basis[i][0],basis[i][1],basis[i][2]);
  fprintf(outfile,"%%endblock LATTICE_CART\n\n");

  fprintf(outfile,"%%block POSITIONS_FRAC\n");
  for(i=0;i<natoms;i++)
    fprintf(outfile,"%3s %.9f %.9f %.9f\n",
                     atno2sym(atoms[i].atno),atoms[i].frac[0],
                     atoms[i].frac[1],atoms[i].frac[2]);
  fprintf(outfile,"%%endblock POSITIONS_FRAC\n");
  cell_write_kpts(outfile);
}


void cell_write_abc(FILE* outfile){
  int i;
  double abc[6];

  fprintf(outfile,"%%block LATTICE_ABC\nang\n");
  cart2abc(basis,abc,1);
  fprintf(outfile,"%f %f %f\n%f %f %f\n",
                  abc[0],abc[1],abc[2],abc[3],abc[4],abc[5]);
  fprintf(outfile,"%%endblock LATTICE_ABC\n\n");

  fprintf(outfile,"%%block POSITIONS_FRAC\n");
  for(i=0;i<natoms;i++)
    fprintf(outfile,"%3s %.9f %.9f %.9f\n",
                     atno2sym(atoms[i].atno),atoms[i].frac[0],
                     atoms[i].frac[1],atoms[i].frac[2]);
  fprintf(outfile,"%%endblock POSITIONS_FRAC\n");
  cell_write_kpts(outfile);
}

void cell_write_abs(FILE* outfile){
  int i;

  fprintf(outfile,"%%block LATTICE_CART\nang\n");

  for(i=0;i<3;i++)
    fprintf(outfile,"%f %f %f\n",
                  basis[i][0],basis[i][1],basis[i][2]);
  fprintf(outfile,"%%endblock LATTICE_CART\n\n");

  fprintf(outfile,"%%block POSITIONS_ABS\n");
  for(i=0;i<natoms;i++)
    fprintf(outfile,"%3s %.9f %.9f %.9f\n",
                     atno2sym(atoms[i].atno),atoms[i].abs[0],
                     atoms[i].abs[1],atoms[i].abs[2]);
  fprintf(outfile,"%%endblock POSITIONS_ABS\n");
  cell_write_kpts(outfile);
}


void cell_write_abc_abs(FILE* outfile){
  int i;
  double abc[6];

  fprintf(outfile,"%%block LATTICE_ABC\nang\n");
  cart2abc(basis,abc,1);
  fprintf(outfile,"%f %f %f\n%f %f %f\n",
                  abc[0],abc[1],abc[2],abc[3],abc[4],abc[5]);
  fprintf(outfile,"%%endblock LATTICE_ABC\n\n");

  fprintf(outfile,"%%block POSITIONS_ABS\n");
  for(i=0;i<natoms;i++)
    fprintf(outfile,"%3s %.9f %.9f %.9f\n",
                     atno2sym(atoms[i].atno),atoms[i].abs[0],
                     atoms[i].abs[1],atoms[i].abs[2]);
  fprintf(outfile,"%%endblock POSITIONS_ABS\n");
  cell_write_kpts(outfile);
}

void cell_write_kpts(FILE* outfile){
  int i,j;
  double v[3];

  if ((symops)&&(nsym>1)){
    fprintf(outfile,"\n%%block SYMMETRY_OPS");
    for(i=0;i<nsym;i++){
      fprintf(outfile,"\n");
      if (debug){
	fprintf(outfile,"! ");
	ident_sym(&symops[i],outfile);
      }
      for(j=0;j<3;j++)
	fprintf(outfile,"%13.9f %13.9f %13.9f\n",symops[i].mat[j][0],
		symops[i].mat[j][1],symops[i].mat[j][2]);
      /* Though the matrix is in cartesian co-ords, the associated
       * translation is expected in fractional co-ords... */
      if (symops[i].tr){
        for(j=0;j<3;j++)
          v[j]=symops[i].tr[0]*recip[j][0]+symops[i].tr[1]*recip[j][1]+
            symops[i].tr[2]*recip[j][2];
        fprintf(outfile,"%13.9f %13.9f %13.9f\n",v[0],v[1],v[2]);
      }
      else
	fprintf(outfile,"0.0 0.0 0.0\n");
    }
    fprintf(outfile,"%%endblock SYMMETRY_OPS\n");
  }

  if ((mp)&&(mp->grid[0]>0)){
    fprintf(outfile,"\nKPOINT_MP_GRID %d %d %d\n",mp->grid[0],
            mp->grid[1],mp->grid[2]);

    fprintf(outfile,"KPOINT_MP_OFFSET %.9f %.9f %.9f\n",mp->disp[0],
            mp->disp[1],mp->disp[2]);
  }
  else if (kpts) {
    fprintf(outfile,"\n%%block KPOINT_LIST\n");
    for(i=0;i<*nkpts;i++)
      fprintf(outfile,"%13.9f %13.9f %13.9f     %12.9f\n",kpts[i].k[0],
              kpts[i].k[1],kpts[i].k[2],kpts[i].wt);
    fprintf(outfile,"%%endblock KPOINT_LIST\n");
  }

}
