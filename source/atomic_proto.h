/* atomicdata.c */
int get_atomic_data(char masterfile[]);
/* atomicdata_sub.c */
int atomicdata2file(void);
int index_lines(void);
int index_phot_top(void);
int index_inner_cross(void);
void indexx(int n, float arrin[], int indx[]);
int limit_lines(double freqmin, double freqmax);
int check_xsections(void);
double q21(struct lines *line_ptr, double t);
double q12(struct lines *line_ptr, double t);
double a21(struct lines *line_ptr);
double upsilon(int n_coll, double u0);
void skiplines(FILE *fptr, int nskip);
/* atomicdata_init.c */
TopPhotPtr alloc_phot_top_record(int n, int np);
TopPhotPtr alloc_inner_cross_record(int n, int np);
AugerPtr alloc_auger_macro_record(int n);
int free_xsection_data(void);
int init_atomic_data(void);
