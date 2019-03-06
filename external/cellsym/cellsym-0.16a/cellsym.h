#define CS_VER "0.16a"

/* Global variables for system description */

#define BOHR 0.529177249

/* constants for the file formats which can be written */
#define PDB 3
#define CELL 4
#define CELL_ABC 5
#define CELL_ABS 6
#define CELL_ABC_ABS 7
#define CNULL 13
#define SHELX 14

/* flags for reading and output */

#define PDBN 512
#define SHELX_AIRSS 2048

/* And for laziness when reading .cell or .pdb file */
#define MAX_ATOMS 2000


struct atom
   {unsigned int atno; double abs[3]; double frac[3]; double force[3];};
struct grid {char *name; int size[3]; double *data; struct grid *next;};
struct kpt {double k[3]; double wt;};
struct mp_grid {int grid[3]; double disp[3];};
struct vector {double v[3]; double mod2;};
struct sym_op {double mat[3][3]; double *tr;};

void error_exit(char* msg);
void real2rec(void);
void addfrac(void);
void addabs(void);
void abc2cart(double *abc, double basis[3][3]);
void basis2abc(double basis[3][3], double abc[6]);
void cart2abc(double basis[3][3], double *abc, int fix);

unsigned int atsym2no(char* sym);
char* atno2sym(unsigned no);

#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#endif

#define aeq(a,b) (fabs((a)-(b))<((fabs(a)+0.5)*tol))

#ifndef CELLSYM_MAIN
extern double (*basis)[3],recip[3][3]; /* Basis sets */
extern double cell_vol;
extern int natoms,nsym,debug,flags;
extern struct atom *atoms;
extern double tol;
extern int *nkpts;
extern struct kpt *kpts;
extern struct mp_grid *mp;
extern int preserve_c;
extern struct sym_op *symops;
extern char *title;
#endif
