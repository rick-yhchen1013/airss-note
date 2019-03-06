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

/* Changed 4/07 to start atom numbering at one, not zero
 * and to put atomic symbols in the correct place */

/* Excerpted from PDB file format (www.pdb.org):
 *
 * Columns    Thing
 *
 *  1-6       "ATOM  "
 *  7-11      Atom serial no (int)
 * 13-16      Atom name (lable)
 *  ...
 * 31-38      x co-ord, as Fortran's real(8.3), Angstroms
 * 39-46      y co-ord, ditto
 * 47-54      z co-ord, ditto
 *  ...
 * 77-78      Element symbol, right justified
 *
 *  1-6       "CRYST1"
 *  7-15      a, as Fortran's real(9.3), Angstroms
 * 16-24      b, ditto
 * 25-33      c, ditto
 * 34-40      alpha, as Fortran's real(7.2), degrees
 * 41-47      beta, ditto
 * 48-54      gamma, ditto
 * 56-66      space group, left justified string
 * 67-70      number of polmeric chains per cell (z)
 *
 * Set z to one, and space group to "P 1" if irrelevant / undetermined
 */

#include<stdio.h>
#include<string.h>
#include<stdlib.h>

#include "cellsym.h"

#define MAX_ELS 103

void pdb_write(FILE* outfile){
  int i,*n_in_el=NULL;
  double abc[6];

  if (flags&PDBN){
    n_in_el=calloc(MAX_ELS+1,sizeof(int));
    if (!n_in_el) error_exit("Calloc error in pdbn_write");
  }

  fprintf(outfile,"REMARK written by check2xsf\n");

  if (title)
    fprintf(outfile,"TITLE %s\n",title);

  cart2abc(basis,abc,1);
  fprintf(outfile,"CRYST1 %8.3f %8.3f %8.3f %6.2f %6.2f %6.2f P 1"
                  "           1\n",
                  abc[0],abc[1],abc[2],abc[3],abc[4],abc[5]);

/* Ditch chain identifier (X)? */

  for(i=0;i<natoms;i++)
    if (flags&PDBN){
      fprintf(outfile,"ATOM   %4d ",i+1);
      if (++n_in_el[min(atoms[i].atno,MAX_ELS)]<=9)
          fprintf(outfile,"%3s%d",atno2sym(atoms[i].atno),
                  n_in_el[min(atoms[i].atno,MAX_ELS)]);
      else if(n_in_el[min(atoms[i].atno,MAX_ELS)]<=99)
          fprintf(outfile,"%2s%d",atno2sym(atoms[i].atno),
                  n_in_el[min(atoms[i].atno,MAX_ELS)]);
      else if((n_in_el[min(atoms[i].atno,MAX_ELS)]<=999)&&
              (strlen(atno2sym(atoms[i].atno))==1))
          fprintf(outfile,"%1s%d",atno2sym(atoms[i].atno),
                  n_in_el[min(atoms[i].atno,MAX_ELS)]);
      else fprintf(outfile,"%3s*",atno2sym(atoms[i].atno));
      fprintf(outfile,"   X     1     %7.3f %7.3f %7.3f"
                   "                      %2s\n",
                     atoms[i].abs[0],
                     atoms[i].abs[1],atoms[i].abs[2],atno2sym(atoms[i].atno));
    }
    else fprintf(outfile,"ATOM   %4d %4s   X     1     %7.3f %7.3f %7.3f"
                   "                      %2s\n",
                     i+1,atno2sym(atoms[i].atno),atoms[i].abs[0],
                     atoms[i].abs[1],atoms[i].abs[2],atno2sym(atoms[i].atno));

  fprintf(outfile,"END\n");
  if (flags&PDBN) free(n_in_el);
}
