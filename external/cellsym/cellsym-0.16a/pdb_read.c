/* A very simplistic .pdb reader */

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
#include<stdlib.h> /* malloc */
#include<string.h>
#include<errno.h>

#include "cellsym.h"

#define LINE_SIZE 100

static int strnncpy(char *dest, char *src, int n);

void pdb_read(FILE* infile){
  double abc[6],*dptr;
  int have_basis=0,i,j;
  char buffer[LINE_SIZE+1],buff2[LINE_SIZE],*cptr;

  natoms=0;

  if (!(atoms=malloc(MAX_ATOMS*sizeof(struct atom))))
       error_exit("Malloc error in pdb_read");
  if (!(basis=malloc(72))) error_exit("Malloc error in pdb_read for basis");

  while(1){
   for(i=0;i<LINE_SIZE-1;i++) buffer[i]=0;
    if (!fgets(buffer,LINE_SIZE,infile)) break;
/* First six characters are record name */
    strnncpy(buff2,buffer,6);

    if (!strcasecmp(buff2,"REMARK")) continue;
    if (!strcasecmp(buff2,"TER")) continue;
    if (!strcasecmp(buff2,"END")) break;
    if (!strcasecmp(buff2,"MTRIX1")) continue; /* MTRXn are symmetry */
    if (!strcasecmp(buff2,"MTRIX2")) continue; /* operations which can */
    if (!strcasecmp(buff2,"MTRIX3")) continue; /* be safely ignored */


    if (!strcasecmp(buff2,"CRYST1")){
      sscanf(buffer+7,"%lf %lf %lf %lf %lf %lf",abc,abc+1,abc+2,
                                              abc+3,abc+4,abc+5);
/* By convention, a=b=c=1, alpha=beta=gamma=90 is a dummy entry here... */
      if ((abc[0]==1)&&(abc[1]==1)&&(abc[2]==1)&&
          (abc[3]==90)&&(abc[4]==90)&&(abc[5]==90)) continue;
      abc2cart(abc,basis);
      have_basis=1;
      continue;
    }

    if ((strcasecmp(buff2,"ATOM")==0)||(strcasecmp(buff2,"HETATM")==0)){
       strnncpy(buff2,buffer+30,24); /* grab co-ords section */
       dptr=atoms[natoms].abs;
       if (sscanf(buff2," %lf %lf %lf",dptr,dptr+1,dptr+2)!=3){
          fprintf(stderr,"Error parsing line\n%s\n",buffer);
          fprintf(stderr,"buff2 was '%s'\n",buff2);
          exit(1);
       }
/* The atomic symbol ought to be in columns 77 to 78 */
       strnncpy(buff2,buffer+76,2);
       if (buff2[0])
         atoms[natoms].atno=atsym2no(buff2);
/* But it might be randomly-justified in cols 13-16 */
       else {
         strnncpy(buff2,buffer+12,4);
         cptr=buff2;
         while(*cptr==' ') cptr++;
         cptr[2]=0;
         if ((cptr[1]>'0')&&(cptr[1]<'9')) cptr[1]=0;
         atoms[natoms].atno=atsym2no(cptr);
       }
       natoms++;
       if (natoms>MAX_ATOMS){
         fprintf(stderr,"Warning: pdb file truncated at %d atoms\n",natoms);
         break;
       }
       continue;
    }

    if (debug)
      fprintf(stderr,"Warning, ignoring PDB line entitled '%s'\n",buff2);
  }

  if(have_basis==0){
    if (debug) fprintf(stderr,"Warning, no unit cell in pdb file\n"
                               "Creating dummy 10A box\n");
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
        basis[i][j]=0;
    for(i=0;i<3;i++) basis[i][i]=10;
  }

  if (debug>1) fprintf(stderr,"%d atoms read\n",natoms);

  real2rec();
  addfrac();

}


/* Function to copy n characters, to null-terminate, and to kill trailing
 * whitespace. Destination must be n+1 characters long, unless source was
 * null terminated */

static int strnncpy(char *dest, char *src, int n){
  int i;
  for (i=0;src[i]&&i<n;i++) dest[i]=src[i];
  dest[i--]=0;
  while((i>=0)&&((dest[i]==' ')||(dest[i]=='\n')||(dest[i]=='\r'))) dest[i--]=0;
  return(i+2);
}
