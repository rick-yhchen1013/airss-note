/* This code glues together some routines from check2xsf with libspg */


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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>  /* isatty */
#include <math.h> /* fabs */
#include <string.h>

#define CELLSYM_MAIN
#include "cellsym.h"
#include <spglib/spglib.h>
// #include <spglib.h>

void help(void);
struct sym_op *sym_frac2abs(int spg_rot[][3][3],double spg_tr[][3],int nsym);

void pdb_write(FILE* outfile);
void cell_write(FILE* outfile);
void cell_write_abc(FILE* outfile);
void cell_write_abs(FILE* outfile);
void cell_write_abc_abs(FILE* outfile);
void shelx_write(FILE* outfile);
void ident_sym(struct sym_op *s, FILE* outfile);

void cell_read(FILE* infile);
void pdb_read(FILE* infile);
void shelx_read(FILE* infile);

/* Global variables for system description */

double (*basis)[3],recip[3][3]; /* Basis sets */
double cell_vol;
int natoms,nsym;
int debug,molecule,format,forces,flags,preserve_c;
struct atom *atoms;
int *nkpts,*symgen;
struct kpt *kpts;
struct mp_grid *mp;
struct sym_op *symops;
double *symtol;
int lhs_fudge;
double tol=1e-4;
char *title;

int main(int argc, char **argv)
{
  int i,j,opt=1,result,primitive,refine,international,schoen;
  int point, symmetry,list,natoms2;
  int dummy[3][3];
  int results,resulti,resultp;
  char *optp,*file1=NULL,*file2=NULL;
  double tolmin,tolmax;
  FILE *infile,*outfile;

  /* SPG variables */
  double spg_latt[3][3];
  double (*spg_pos)[3],spg_symprec;
  int *spg_type;
  char spg_sym[22];
  SpglibDataset *spg;

/* Initialise all pointers to NULL, etc. */

  atoms=NULL;
  basis=NULL;
  nkpts=NULL;
  kpts=NULL;
  mp=NULL;
  title=NULL;
  spg=NULL;
  debug=0;
  preserve_c=0;
  natoms=0;
  nsym=0;
  flags=0;
  tolmin=tolmax=1e-4;
  for(i=0;i<3;i++){
    for(j=0;j<3;j++) recip[i][j]=0.0;
  }

  primitive=refine=international=schoen=point=symmetry=list=0;
  format=0;
  opt=1;
  optp=argv[opt];

  while (opt<argc){
    switch(*optp){
    case 0:
      opt++;
      optp=argv[opt];
      break;
    case '-':
      if(*(optp+1)=='-'){ /* Gnu-style "--" option */
	if (!strcmp(optp,"--pdb")) format=PDB;
        else if (!strcmp(optp,"--pdbn")) {format=PDB; flags|=PDBN;}
        else if (!strcmp(optp,"--cell")) format=CELL;
        else if (!strcmp(optp,"--cell_abc")) format=CELL_ABC;
        else if (!strcmp(optp,"--cell_abs")) format=CELL_ABS;
        else if (!strcmp(optp,"--cell_abc_abs")) format=CELL_ABC_ABS;
        else if (!strcmp(optp,"--shelx")) format=SHELX;
        else if (!strcmp(optp,"--airss")) {format=SHELX; flags|=SHELX_AIRSS;}
        else if (!strcmp(optp,"--null")) format=CNULL;
	else if (!strcmp(optp,"--primitive")) primitive=1;
	else if (!strcmp(optp,"--refine")) refine=1;
	else if (!strcmp(optp,"--int")) international=1;
	else if (!strcmp(optp,"--schoen")) schoen=1;
	else if (!strcmp(optp,"--symmetry")) symmetry=1;
	else if (!strcmp(optp,"--point")) point=1;
	else if (!strcmp(optp,"--list")) list=1;
        else if (!strcmp(optp,"--help")) help();
       else {
          fprintf(stderr,"Invalid option %s.\n%s -h for usage.\n",
                   optp,argv[0]);
          exit(1);
        }
        opt++;
        optp=argv[opt]; 
        break;
      }
      while(*(++optp)){
        switch(*optp){
        case 'v':
          debug++;
          break;
        case 'h':
          help();
          break;
        case 'V':
          printf("cell2sym version " CS_VER "\n");
#ifdef SPGLIBVER
          printf("Linked with spglib version %d.%d.%d\n",
                 spg_get_major_version(),
                 spg_get_minor_version(),
                 spg_get_micro_version());
#endif
          exit(0);
          break;
        case 'e':
          if(*(optp+1)=='='){
	    i=sscanf(optp+2,"%lf-%lf",&tolmin,&tolmax);
	    if (i==1) tolmax=tolmin;
	    else if (i!=2) error_exit("malformed option -e=");
          }
          else error_exit("malformed option -e");
          while(*((++optp)+1));
          break;
        default:
          fprintf(stderr,"Invalid option %c.\n%s -h for usage.\n",
                   *optp,argv[0]);
           exit(1);
         }
       }
       opt++;
       optp=argv[opt];
       break;
    default:
      if (!file1) file1=argv[opt];
      else if (!file2) file2=argv[opt];
      else{
	fprintf(stderr,"Unexpected argument %s\n%s -h for usage.\n",
		argv[opt],argv[0]);
	exit(1);
      }
      opt++;
      optp=argv[opt];
      break;
    }
  }

  if (debug>1) fprintf(stderr,"cellsym version " CS_VER 
#ifdef SPGLIBVER
          " linked with spglib version %d.%d.%d\n",
                 spg_get_major_version(),
                 spg_get_minor_version(),
                 spg_get_micro_version()
#else
         "\n"
#endif
  );

  if (!file1) error_exit("no input file specified.");

  if (format==0){
    if((primitive==0)&&(refine==0)&&(symmetry==0)&&(file2==NULL)) format=CNULL;
    else format=CELL;
  }

  if (list) symmetry=1;

  if ((!file2)&&(isatty(fileno(stdout)))&&(format!=CNULL))
    error_exit("refusing to output to a terminal");

  infile=fopen(file1,"rb");
  if(!infile){
    fprintf(stderr,"Error, unable to open %s for reading.\n",file1);
    exit(1);
  }

  if (file2){
    outfile=fopen(file2,"wb");
    if(!outfile){
      fprintf(stderr,"Error, unable to open %s for writing.\n",file2);
      exit(1);
    }
  }
  else
    outfile=stdout;

  i=strlen(file1);
  if((i>4)&&(!strcmp(file1+i-4,".pdb")))
    pdb_read(infile);
  if((i>4)&&(!strcmp(file1+i-4,".res")))
    shelx_read(infile);
  else
    cell_read(infile);
  fclose(infile);

  if (cell_vol<0) cell_vol=fabs(cell_vol);

  /* Now fill in the SPG-style variables */

  if (refine){
    spg_pos=malloc(4*3*natoms*sizeof(double));
    spg_type=malloc(4*natoms*sizeof(int));
  }
  else{
    spg_pos=malloc(3*natoms*sizeof(double));
    spg_type=malloc(natoms*sizeof(int));
  }

  results=resulti=resultp=-1;
  natoms2=natoms;
  tol=tolmin;
  while(tol<=tolmax){

    spg_symprec=tol;
    tol*=2;

    for(i=0;i<natoms;i++){
      for(j=0;j<3;j++)
	spg_pos[i][j]=atoms[i].frac[j];
      spg_type[i]=atoms[i].atno;
    }


    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	spg_latt[i][j]=basis[j][i];

    if (primitive){
      natoms2=spg_find_primitive(spg_latt,spg_pos,spg_type,natoms,spg_symprec);
    }

    if (refine){
      natoms2=spg_refine_cell(spg_latt,spg_pos,spg_type,natoms,spg_symprec);
    }

    if (schoen){
      results=spg_get_schoenflies(spg_sym,spg_latt,spg_pos,spg_type,natoms,
				 spg_symprec);

      if (!result)
	fprintf(stderr,"Attempt to find Schoenflies symmetry failed\n");
      else if (result!=results){
	fprintf(stderr,"Tol=%g Schoenflies symmetry is %s\n",
		spg_symprec,spg_sym);
	results=result;
      }
    }
  
    if (symmetry||point||international){
      spg=spg_get_dataset(spg_latt,spg_pos,spg_type,natoms,spg_symprec);
      if (!spg){
        fprintf(stderr,"Attempt to find symmetry failed\n");
        exit(1);
      }
      if (international){
	result=spg->spacegroup_number;
	if (!result)
	  fprintf(stderr,"Attempt to find international symmetry failed\n");
	else if (result!=resulti){
	  fprintf(stderr,"Tol=%g International symmetry is %s\n",
		  spg_symprec,spg->international_symbol);
	  resulti=result;
	}
      }
      if (symmetry||point){
	result=spg->n_operations;
	if (result!=nsym){
	  nsym=result;
	  if (debug) fprintf(stderr,"Tol=%g %d symmetry operations found.\n",
			     spg_symprec,nsym);
	}
	if (point){
	  result=spg_get_pointgroup(spg_sym,dummy,spg->rotations,nsym);
	  if (!result)
	    fprintf(stderr,"Attempt to find point group failed\n");
	  else if (result!=resultp){
	    fprintf(stderr,"Tol=%g Point group is %s\n",spg_symprec,spg_sym);
	    resultp=result;
	  }
	}
      }
    }
  }
    
  natoms=natoms2;

  /* Convert back to c2x-style variables */
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      basis[j][i]=spg_latt[i][j];
  real2rec();

  free(atoms);
  atoms=malloc(natoms*sizeof(struct atom));
  for(i=0;i<natoms;i++){
    for(j=0;j<3;j++)
      atoms[i].frac[j]=spg_pos[i][j];
    atoms[i].atno=spg_type[i];
  }
  addabs();

  /* SPG has sym-ops in fractional co-ords, Castep in cartesians */
  if (nsym) {
    symops=sym_frac2abs(spg->rotations,spg->translations,nsym);
    if (list)
      for(i=0;i<nsym;i++) ident_sym(&symops[i],stderr);
  }

  switch(format){
    case PDB:
      pdb_write(outfile);
      break;
    case CELL:
      cell_write(outfile);
      break;
    case CELL_ABC:
      cell_write_abc(outfile);
      break;
    case CELL_ABS:
      cell_write_abs(outfile);
      break;
    case CELL_ABC_ABS:
      cell_write_abc_abs(outfile);
      break;
    case SHELX:
      shelx_write(outfile);
      break;
    case CNULL:
      break;

    default:
      fprintf(stderr,"This cannot happen. Sorry\n");
  }

  return(0);

}

struct sym_op *sym_frac2abs(int spg_rot[][3][3],double spg_tr[][3],int nsym){
  int i,j,k,kk;
  struct sym_op *s;
  double mat[3][3];

  s=malloc(nsym*sizeof(struct sym_op));
  if(!s) error_exit("Malloc error in sym_frac2abs");


  for(i=0;i<nsym;i++){
    if (debug>2){
      fprintf(stderr,"Sym op on entry:\n");
      for(j=0;j<3;j++)
	fprintf(stderr,"%2d %2d %2d\n",spg_rot[i][j][0],
	      spg_rot[i][j][1],spg_rot[i][j][2]);
      fprintf(stderr,"%13.9f %13.9f %13.9f\n",spg_tr[i][0],spg_tr[i][1],
	      spg_tr[i][2]);
    }
    /* Translations are "easy" */
    if ((spg_tr[i][0]!=0)||(spg_tr[i][1]!=0)||(spg_tr[i][2]!=0)){
      s[i].tr=malloc(3*sizeof(double));
      if(!s[i].tr) error_exit("Malloc error for tr in sym_frac2abs");
      for(j=0;j<3;j++) s[i].tr[j]=spg_tr[i][0]*basis[0][j]+
			 spg_tr[i][1]*basis[1][j]+spg_tr[i][2]*basis[2][j];
    }
    else s[i].tr=NULL;
    /* Matrices are harder. Permute indices until answer looks right? */
    for(j=0;j<3;j++)
      for(k=0;k<3;k++)
	mat[j][k]=s[i].mat[j][k]=0.0;
    for(j=0;j<3;j++)
      for(k=0;k<3;k++)
	for(kk=0;kk<3;kk++)
	  mat[j][k]+=spg_rot[i][j][kk]*recip[kk][k];
    for(j=0;j<3;j++)
      for(k=0;k<3;k++)
	for(kk=0;kk<3;kk++)
	  s[i].mat[j][k]+=basis[kk][j]*mat[kk][k];
    if (debug>2){
      fprintf(stderr,"Sym op on exit:\n");
      for(j=0;j<3;j++)
	fprintf(stderr,"%13.9f %13.9f %13.9f\n",s[i].mat[j][0],
                s[i].mat[j][1],s[i].mat[j][2]);
      if (s[i].tr)
	fprintf(stderr,"%13.9f %13.9f %13.9f\n",s[i].tr[0],
                s[i].tr[1],s[i].tr[2]);
      else
	fprintf(stderr,"0.0 0.0 0.0\n");
    }
  }

  return s;
}

void help(void){
  printf("Usage: cellsym [-vV] [-e=eps] [--FORMAT]"
	 " --OPERATION infile [outfile]\n\n"
	 "-e=eps       set tolerance for finding symmetries. Default=%g\n"
	 "-v           be verbose (may be repeated)\n"
	 "-V           print version number and exit\n\n",tol);
  printf("OPERATION is one of:\n"
	 "   primitive      call spg_find_primitive()\n"
	 "   refine         call spg_refine_cell()\n"
	 "   int            call spg_get_dataset()"
	 " and report international symbol\n"
	 "   schoen         call spg_get_schoenflies()\n"
	 "   symmetry       call spg_get_dataset() and keep symmetry ops\n"
	 "   list           call spg_get_dataset() and list symmetry ops\n"
	 "   point          call spg_get_dataset() followed by "
	 "spg_get_pointgroup()\n\n");
  printf("FORMAT is one of: cell      CASTEP .cell, cartesian and fractional\n"
         "                  cell_abc                abc and fractional\n"
         "                  cell_abs                cartesian and absolute\n"
         "                  cell_abc_abs            abc and absolute\n"
         "                  pdb       PDB\n"
         "                  pdbn      PDB with atoms numbered\n"
         "                  shelx     SHELX97\n"
         "                  airss     not SHELX97\n\n");
  printf("Multiple operations may be specified. The order in which they will\n"
	 "be applied is the order in which they are listed above, not the\n"
	 "order on the command line.\n\n");
  printf("cellsym version " CS_VER
	 " (c) 2012-2016 MJR, symmetry information from libspg.\n\n");
  exit(0);
}

void error_exit(char *msg){
  fprintf(stderr,"Aborting: %s\n",msg);
  exit(1);
}
