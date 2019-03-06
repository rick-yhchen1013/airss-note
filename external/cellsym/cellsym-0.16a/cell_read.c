/* Read some useful data from a CASTEP .cell file
 * 
 * Should read lattice_cart, lattice_abc, positions_frac and positions_abs
 * blocks. Should skip blanks lines and comments, and should cope with
 * UNIX, DOS or Mac line-endings.
 */


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

static int cellreadline(char *buffer, int len, FILE* infile);

static int line_count;
static char *ptr;

void cell_read(FILE* infile){
  int units,i,j,frac;
  double lat_abc[6],*dptr;
  char sym[4],*ptr2;
  static char buffer[LINE_SIZE+1];

  line_count=0;
  units=0;
  natoms=0;
  frac=0;

  if (debug>2) fprintf(stderr,"Cell read called\n");

  if (!(basis=malloc(72))) error_exit("Malloc error in cellread for basis");

  while(cellreadline(buffer,LINE_SIZE,infile)){
    ptr=buffer;
/* What remains ought to be a block or a keyword. We care for blocks only */
    if (strncasecmp(ptr,"%block",6)) continue;
    ptr+=6;
/* Eat leading spaces */
    while(*ptr==' ') ptr++;
/* Kill trailing spaces and things after ! */
    ptr2=ptr;
    while((*ptr2)&&(*ptr2!=' ')&&(*ptr2!='!')) ptr2++;
    *ptr2=0; /* Either it was a null, or it should be one... */

    if (debug>2) fprintf(stderr,"Found a %s block\n",ptr);

    if (!strcasecmp(ptr,"lattice_cart")){
      cellreadline(buffer,LINE_SIZE,infile);
      ptr=buffer;
      if (!strncasecmp(ptr,"ang",3)){
        units=0;
        cellreadline(buffer,LINE_SIZE,infile);
      }else if(!strncasecmp(ptr,"bohr",4)){
        units=1;
        cellreadline(buffer,LINE_SIZE,infile);
      }
      for(i=0;i<3;i++){
        if(sscanf(buffer,"%lf %lf %lf",basis[i],basis[i]+1,basis[i]+2)!=3){
          fprintf(stderr,"Parse error in cell_read for line %d\n",line_count);
           exit(1);
        }
        cellreadline(buffer,LINE_SIZE,infile);
      }
      if (units==1)
        for(i=0;i<3;i++)
          for(j=0;j<3;j++)
            basis[i][j]*=BOHR;

    }else if(!strcasecmp(ptr,"lattice_abc")){
      cellreadline(buffer,LINE_SIZE,infile);
      ptr=buffer;
      if (!strncasecmp(ptr,"ang",3)){
        units=0;
        cellreadline(buffer,LINE_SIZE,infile);
      }else if(!strncasecmp(ptr,"bohr",4)){
        units=1;
        cellreadline(buffer,LINE_SIZE,infile);
      }
      for(i=0;i<2;i++){
        if(sscanf(buffer,"%lf %lf %lf",lat_abc+3*i,lat_abc+3*i+1,
                                    lat_abc+3*i+2)!=3){
          fprintf(stderr,"Parse error in cell_read for line %d\n",line_count);
          if (debug) fprintf(stderr,"%s\n",buffer);
          exit(1);
        }
        cellreadline(buffer,LINE_SIZE,infile);
      }
      if (units==1)
        for(i=0;i<3;i++)
          *(lat_abc+i)*=BOHR;
      abc2cart(lat_abc,basis);

    }else if((!strcasecmp(ptr,"positions_frac"))||
            (!strcasecmp(ptr,"positions_abs"))){
      frac=1;
      if(!strcasecmp(ptr,"positions_abs")) frac=0;
      cellreadline(buffer,LINE_SIZE,infile);
      ptr=buffer;
      units=0;
      if(!frac){
        if (!strncasecmp(ptr,"ang",3)){
          cellreadline(buffer,LINE_SIZE,infile);
        }else if(!strncasecmp(ptr,"bohr",4)){
          units=1;
          cellreadline(buffer,LINE_SIZE,infile);
        }
      }
      if (!(atoms=malloc(MAX_ATOMS*sizeof(struct atom))))
          error_exit("Malloc error in cell_read");
      for(i=0;i<MAX_ATOMS;i++){
        if (frac) dptr=atoms[i].frac;
        else dptr=atoms[i].abs;
        if(sscanf(buffer,"%d %lf %lf %lf",&atoms[i].atno,dptr,
                                      dptr+1,dptr+2)==4);
        else if(sscanf(buffer,"%3s %lf %lf %lf",sym,dptr,
                                      dptr+1,dptr+2)==4)
          atoms[i].atno=atsym2no(sym);
        else if(!strncasecmp(ptr,"%endblock",9)) break;
        else{
          fprintf(stderr,"*Parse error in cell_read for line %d\n",line_count);
          if (debug) fprintf(stderr,"%s\n",ptr); 
          exit(1);
        }
        cellreadline(buffer,LINE_SIZE,infile);
      }
      natoms=i;
      if (debug>2) fprintf(stderr,"%d atoms read\n",natoms);
      if (units==1)
        for(i=0;i<natoms;i++)
          for(j=0;j<3;j++)
            atoms[i].abs[j]*=BOHR;
    }else{
      while(cellreadline(buffer,LINE_SIZE,infile)){
        ptr=buffer;
        if (!strncasecmp(ptr,"%endblock",9)) break;
      }
    }
  }

  real2rec();
  if (frac) addabs();
  else addfrac();

}

int cellreadline(char *buffer, int len, FILE* infile){
  int off;
  char *ptr,*success,*p2;

  while((success=fgets(buffer,len,infile))){ /* fgets() always
                                                    null terminates,
                                                    gcc likes extra brackets */
    line_count++;

/* Kill trailing spaces and newlines / carriage returns */
    ptr=buffer+strlen(buffer)-1;
    while((ptr>=buffer)&&((*ptr==' ')||(*ptr=='\n')||(*ptr=='\r'))) ptr--;
    *(ptr+1)=0;

/* Eat leading spaces */
    ptr=buffer;
    while(*ptr==' ') ptr++;
/* Parse a very special comment */
    if (!strncasecmp(ptr,"#titl ",6)){ /* We have a shelx-style title */
      ptr+=6;
      p2=ptr;
      while((*p2)&&(*p2!='\n')&(*p2!='\r')) p2++;
      title=malloc(p2-ptr+1);
      if(!title) error_exit("Malloc error in cell_read()");
      strncpy(title,ptr,p2-ptr);
      title[p2-ptr]=0;
      continue;
    }
/* Skip comments and blank lines */
    if (!strncasecmp(ptr,"comment",7)) continue;
    if ((*ptr=='#')||(*ptr=='!')||(*ptr==';')||(*ptr==0)) continue;
    break;
  }

  if (!success) return(0);

  off=ptr-buffer;
  if (off){
    while(*ptr) {*(ptr-off)=*ptr;ptr++;}
    *(ptr-off)=0;
  }

  return (1);
}
