/* A very simplistic SHELX reader */


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


/* Does not handle continuation lines.
 *
 * Parses CELL and SFAC "commands"
 *
 * Recognises all(?) other commands, so what remains must be an atom
 * label, and is parsed as same.
 */

#include<stdio.h>
#include<stdlib.h> /* malloc */
#include<string.h>
#include<errno.h>
#include<ctype.h>

#include "cellsym.h"

#define LINE_SIZE 100

void shelx_read(FILE* infile){
  double abc[6],dummy;
  int have_basis=0,i,c,is_cmd;
  int spec[102],nspec=1;
  char buffer[LINE_SIZE+1],*cptr,*cptr2;

  char *cmd[]={"TITL","CELL","ZERR","LATT","SFAC","UNIT",
	       "SYMM","DISP","LAUE","REM","MORE","TIME","END","HKLF",
	       "OMIT","SHEL","BASF","TWIN","EXTI","SWAT","HOPE","MERG",
	       "SPEC","RESI","MOVE","ANIS","AFIX","HFIX","FRAG","FEND",
	       "EXYZ","EADP","EQIV","CONN","PART","BIND","DFIX","DANG",
	       "BUMP","SAME","SADI","CHIV","FLAT","DELU","SIMU","DEFS",
	       "ISOR","NCSY","SUMP","L.S.","CGLS","BLOC","DAMP","STIR",
	       "WGHT","FVAR","BOND","CONF","MPLA","RTAB","HTAB","LIST",
	       "ACTA","SIZE","TEMP","WPDB","FMAP","GRID","PLAN","MOLE",0};

  natoms=0;

  if (!(atoms=malloc(MAX_ATOMS*sizeof(struct atom))))
       error_exit("Malloc error in shelx_read");
  if (!(basis=malloc(72))) error_exit("Malloc error in shelx_read for basis");

  while(1){
    for(i=0;i<LINE_SIZE;i++) buffer[i]=0;
    if (!fgets(buffer,LINE_SIZE-1,infile)) break;

    cptr=buffer;

    if (*cptr==0) continue;   /* Skip blank lines */
    if (*cptr=='!') continue; /* Skip comments */
    if (*cptr==' ') continue; /* Lines starting with spaces are comments */

    cptr2=cptr;
    while((*cptr2)&&(*cptr2!=' ')&&(*cptr2!='_')&&(*cptr2!='\n')
	  &&(*cptr2!='\r')&&((cptr2-cptr)<4)) cptr2++;
    *cptr2=0;

    if (!strcasecmp(cptr,"REM")) continue;
    if (!strcasecmp(cptr,"TITL")){
      cptr2++;
      cptr=cptr2;
      while((*cptr2)&&(*cptr2!='\n')&(*cptr2!='\r')) cptr2++;
      title=malloc(cptr2-cptr+1);
      if(!title) error_exit("Malloc error in shelx_read()");
      strncpy(title,cptr,cptr2-cptr);
      title[cptr2-cptr]=0;
      continue;
    }
    if (!strcasecmp(cptr,"CELL")){
      c=sscanf(cptr2+1,"%lf %lf %lf %lf %lf %lf %lf",&dummy,
	       abc,abc+1,abc+2,abc+3,abc+4,abc+5);
      if (c!=7) fprintf(stderr,cptr+4);
      if (c!=7) error_exit("error parsing CELL line");
      have_basis=1;
      continue;
    }
    if (!strcasecmp(cptr,"SFAC")){ /* We now have an unknown number of
				    * atomic species */
      cptr=cptr2+1;
      while(*cptr!=0){
	while(*cptr==' ') cptr++;
	cptr2=cptr;
	while(isalpha(*cptr2)) cptr2++;
	*cptr2=0;
	spec[nspec]=atsym2no(cptr);
	cptr=cptr2+1;
	nspec++;
      }
      if (debug) fprintf(stderr,"%d species in shelx file\n",nspec-1);
      continue;
    }

    /* Do we have another command? */
    is_cmd=0;
    for(i=0;cmd[i];i++){
      if (!strcasecmp(cptr,cmd[i])){
	is_cmd=1;
	break;
      }
    }

    if(!is_cmd){ /* We have an atom */
      c=sscanf(cptr2+1,"%d %lf %lf %lf",&i,atoms[natoms].frac,
	       atoms[natoms].frac+1,atoms[natoms].frac+2);
      if (c!=4){
	fprintf(stderr,"Error parsing following line as atom.\n%s\n",buffer);
	fprintf(stderr,"Blundering on regardless\n");
      }
      else{
	if (i>=nspec) error_exit("Species number too big in shelx_read");
	atoms[natoms].atno=spec[i];
	natoms++;
      }
    }
  }

  if (!have_basis) error_exit("No basis in shelx file");

  abc2cart(abc,basis);
 
  if (debug>1) fprintf(stderr,"%d atoms read\n",natoms);

  real2rec();
  addabs();

}
