/*
 *
 *  ENVIRONMENTAL MODELLING SUITE (EMS)
 *  
 *  File: model/hd-us/include/proto.h
 *  
 *  Description:
 *  Include file for all routine
 *  prototypes
 *  
 *  Copyright:
 *  Copyright (c) 2018. Commonwealth Scientific and Industrial
 *  Research Organisation (CSIRO). ABN 41 687 119 230. All rights
 *  reserved. See the license file for disclaimer and full
 *  use/redistribution conditions.
 *  
 *  $Id: proto.h 6459 2020-02-18 23:42:59Z her127 $
 *
 */

#ifndef _EMSHD_PROTOS_H
#define _EMSHD_PROTOS_H

/*------------------------------------------------------------------*/
/* Version information                                              */
/*------------------------------------------------------------------*/
#define COMPAS_MAJOR_VERSION 1
#define COMPAS_MINOR_VERSION 1
#define COMPAS_PATCH_VERSION 3

/*------------------------------------------------------------------*/
/* Parameter input routines                                         */
/*------------------------------------------------------------------*/
parameters_t *params_alloc(void);
parameters_t *params_read(FILE * fp);
parameters_t *params_read_t(FILE * fp);
parameters_t *auto_params(FILE * fp, int autof);
void auto_params_roam_pre1(FILE * fp, parameters_t *params);
void auto_params_roam_pre2(FILE * fp, parameters_t *params);
void auto_params_roam_pre3(FILE * fp, parameters_t *params);
void auto_params_recom_pre1(FILE * fp, parameters_t *params);
void auto_params_recom_pre2(FILE * fp, parameters_t *params);
void auto_params_roam_post1(FILE * fp, parameters_t *params);
void auto_params_roam_post2(FILE * fp, parameters_t *params);
void auto_params_roam_post3(FILE * fp, parameters_t *params);
void auto_params_roam_post4(FILE * fp, parameters_t *params);
void auto_params_roam_post5(FILE * fp, parameters_t *params);
void auto_params_roam_post6(FILE * fp, parameters_t *params);
void auto_params_recom_post1(FILE * fp, parameters_t *params);
void auto_params_recom_post2(FILE * fp, parameters_t *params);
double get_restart_time(char *filename, char *iunits);
double get_nrt_time(char *filename, char *iunits);
void readparam(FILE * fp);
void read_hdiff(parameters_t *params, FILE *fp, int mode);
int read_blocks(FILE *fp, char *key, int *nb, int **listi, int **listj);
void set_default_param(parameters_t *params);
void params_write(parameters_t *params, dump_data_t *dumpdata);
void cookie_cut(master_t *master, parameters_t *params);
void z0_init(parameters_t *params);
void create_avhrr_list(parameters_t *params);
void create_ghrsst_list(parameters_t *params);
void read_grid(parameters_t *params);
void read_swr(parameters_t *params, FILE *fp, int mode);
void read_hf_bulk(parameters_t *params, FILE *fp);
void read_totals(parameters_t *params, FILE *fp);
void read_exclude_points(parameters_t *params, FILE *fp);
void read_trfilter(parameters_t *params, FILE *fp);
void read_decorr(parameters_t *params, FILE *fp, int mode);
void read_monotone(parameters_t *params, FILE *fp, int mode);
int read_dhw(parameters_t *params, FILE *fp);
void read_trflux(parameters_t *params, FILE *fp);
void tracer_setup(parameters_t *params, FILE *fp);
int numbers_init(parameters_t *params);
char *otime(double dt, char *tag);
void autoset(parameters_t *params, master_t *master, geometry_t *geom);
void autoset_roam(parameters_t *params, master_t *master, geometry_t **window);
void set_sigma_distrib(parameters_t *params, double bmin, double bmax);
void save_sigma_layers(parameters_t *params, master_t *master);
void testtopo(parameters_t *params, master_t *master);
double **set_bathy(parameters_t *params, int cdfid);
void set_bathy_us(parameters_t *params);
double **bathy_read(parameters_t *params, int cdfid, int ti);
double *bathy_read_us(parameters_t *params, int cdfid, int ti);
int decode_tag(char *list, char *tag, char *value);
int get_smoothing(char *list, char *tag);
double get_scaling(char *list, char *tag);
double get_run_length(FILE *fp);
void read_waves(parameters_t *params, FILE *fp, int mode);
void intp_undef(int nce1, int nce2, double **a);
double *value_init(master_t *master, parameters_t *params, char *keyword, char *vname);
void value_init_3d(master_t *master, double *ret, FILE *fp, char *fname, char *vname,
		      char *tag, double fill, char *i_rule);
void value_init_2d(master_t *master, double *ret, FILE *fp, char *fname, char *vname,
		    char *tag, double fill, char *i_rule);
void value_init_sed(master_t *master, double **ret, FILE *fp, char *fname, char *vname,
		    char *tag, double fill, char *i_rule);
int value_init_regions(master_t *master, char *dname, double *tr, int mode);
void trans_write(hd_data_t *hd_data);
char *trtypename(int m, char *buf);
void bathy_compare(master_t *master);
void read_compatible(parameters_t *params, FILE *fp);
void read_means(parameters_t *params, FILE *fp, int mode);
void read_debug(parameters_t *params, FILE *fp);
void read_profile(parameters_t *params, FILE *fp);
void read_explicit_maps(parameters_t *params, FILE *fp);
void get_output_path(parameters_t *params, FILE *fp);
void vel_init(geometry_t *geom, parameters_t *params, master_t *master);
void delaunay_grid(parameters_t *params);
int iswet(int bathy);
int is_wet(parameters_t *params, double bathy);
int is_land(parameters_t *params, double bathy);
int is_outside(parameters_t *params, double bathy);
int is_edge(parameters_t *params, int cc, double *bathy, int **neic);
int has_wet_neighbour(parameters_t *params, int cc, double *bathy, int **neic);
int has_outside_neighbour(parameters_t *params, int cc, double *bathy, int **neic);
int has_neighbour(int c, int npe, int **neic);
void not_implemented(char *text);
void not_included(char *text);
void neighbour_finder(mesh_t *mesh);
void neighbour_finder_a(int ns2, int *npe, int **vedge, int **neic, int **neij);
void neighbour_finder_b(delaunay* d, int ***neic);

int find_neighbour(int c, double **x, double **y, int *npe, int ns2, int *j);
int find_neighbour_l(int c, int ***eloc, int *npe, int ns2, int *j);
void neighbour_none(dump_data_t *dumpdata);

/*------------------------------------------------------------------*/
/* Grid generation                                                  */
/*------------------------------------------------------------------*/
void geog_false_pole_coord(double **x, double **y, double **h1,
                           double **h2, double **a1, double **a2,
                           long int nce1, long int nce2, double x00,
                           double y00, double flon, double flat,
                           double xinc, double yinc);
void geog_dreckon_coord(double **x, double **y, double **h1, double **h2,
                        double **a1, double **a2, long int nce1,
                        long int nce2, double x00, double y00, double rotn,
                        double xinc, double yinc);
void f_pole(double ilat, double ilon, double ang, double *olat,
            double *olon);
int delaunay_dreckon_coord(double *x, double *y,
			   long int nce1, long int nce2,
			   double x00, double y00, 
			   double rotn, double xinc);
int delaunay_rect_coord(double *x, double *y, long int nce1, long int nce2,
			double x00, double y00,	double rotn, double xinc);
void convert_hex_mesh(parameters_t *params, delaunay *d, int mode);
void convert_hex_mesh_orig(parameters_t *params, delaunay *d);
void convert_tri_mesh(parameters_t *params, delaunay *d);
void convert_quad_mesh(parameters_t *params);
delaunay *create_tri_mesh(int np, point *pin, int ns, int *sin, int nh, double *hin, char *code);
void circle_mesh(parameters_t *params, double x00, double y00, double r);
void square_mesh(parameters_t *params, double x00, double y00, double r);
void do_lloyd(parameters_t *params, int np, point *pin, int ns, int *sin, int imax);
void convert_jigsaw_msh(parameters_t *params, char *infile,
			void *jmsh);
void mesh_ofile_name(parameters_t *params, char *key);
delaunay *make_dual(geometry_t *geom, int kin);
int mesh_expand_w(geometry_t *window, int *vec);
int mesh_expand_3d(geometry_t *window, double *u1);

/*------------------------------------------------------------------*/
/* Window subroutines                                               */
/*------------------------------------------------------------------*/
void window_build(geometry_t *geom, parameters_t *params);
void get_local_maps(geometry_t *geom, geometry_t *window, int wn, int *wsa, int wsize,
                    int *ws2, int wsizes2D, int cellf);

void local_map_build(int c, int cc, int wn, int *map, int *nmap,
                     int *rmap, int *wsa, int *ac, int nsaux, int *zoomf);
void local_map_build_c2c(geometry_t * geom, int c, int cc, int wn, int **map, int **nmap,
			 int *wsa, int *ac, int nsaux, int npe);
void local_map_build_z(geometry_t * geom, int c, int cc, int wn, int *map, int *nmap,
		       int *rmap, int *wsa, int *ac, int nsaux);
void local_map_build_o(int c, int cc, int wn, int *map, int *nmap,
                     int *rmap, int *wsa, int *ac, int nsaux, int zoomf);
void get_local_wsc(geometry_t *geom, int *vec, int nvec, int nvec2D, geometry_t **window,
                   int nwindows);
void get_local_wse(geometry_t *geom, int *vec, int nvec, int nvec2D, geometry_t **window,
                   int nwindows, double smag);
void get_local_wsv(geometry_t *geom, int *vec, int nvec, int nvec2D, geometry_t **window,
                   int nwindows);
void get_window_cells(geometry_t *geom, int wn, int *wsa, int *wsize,
		      int *ws2, int wsizeS);
void get_window_cells_h(geometry_t *geom, int wn, int *wsa, int *wsize,
                        int *ws2, int wsizeS);
void get_gl_maps(geometry_t *geom, int nwindows, int **ws2, int *wsizeS);
void set_mask(geometry_t *window);
void window_cells_zoom(geometry_t *geom, parameters_t *params,
			int nwindows, int *zmfe1, int *zmfe2, int **ws2, 
			int *wsizeS);
void window_cells_grouped(geometry_t *geom, int nwindows, int **ws2, int *wsizeS);
void window_cells_region(geometry_t *geom, int nwindows, int **ws2, int *wsizeS, char *fname);
void window_cells_metis(geometry_t *geom, int nwindows, int **ws2, int *wsizeS);
void window_cells_linear_e1(geometry_t *geom, int nwindows, int **ws2,
			    int *wsizeS);
void window_cells_linear_e2(geometry_t *geom, int nwindows, int **ws2,
			    int *wsizeS);
void window_cells_block_e1(geometry_t *geom, int nwindows, int **ws2,
			   int *wsizeS, int os);
void window_cells_block_e2(geometry_t *geom, int nwindows, int **ws2,
			   int *wsizeS, int os);
int reset_map(geometry_t *geom, int cm, int wn, int *wsa, int *al);
void reorder_gl_map(geometry_t *geom, int wn, int *wsa, int wsize, int *mask, int msize);
void wsa_cells(geometry_t *window, int *loc3, int *loc2, int c,
               int *mask, int mode);
void surfbot_build(geometry_t *geom, geometry_t *window, int wn, int *sur, int *bot,
                   int *vec, int nvec, int *vec2D, int nvec2D,
                   int avec2D);
void surfbot_builde(geometry_t *geom, geometry_t *window, int wn, int *sur, int *bot,
                   int *vec, int nvec, int *vec2D, int nvec2D,
                   int avec2D);
geometry_t *window_alloc(void);
void window_init(geometry_t *geom, geometry_t **window);
void win_geom_clear(geometry_t **window, int nwindows);
window_t *win_data_alloc(void);
window_t **win_data_build(master_t *master, geometry_t **window);
window_t *win_data_init(master_t *master, geometry_t *window);
void win_data_clear(window_t *windat);
win_priv_t *win_consts_alloc(void);
win_priv_t **win_consts_init(master_t *master, geometry_t **window);
void win_consts_clear(geometry_t **window, int nwindows);

void get_timesteps(geometry_t **window, window_t **windat,
                   win_priv_t **wincon, int nwindows, master_t *master);
void fill_multidt_aux(master_t *master, geometry_t *window,
                      window_t *windat);
void fill_multidt_aux_2d(geometry_t *window, window_t *windat, double *vel,
                         double *ae);
void timeaux(geometry_t **window, window_t **windat, win_priv_t **wincon,
             int nwindows);
void psl(geometry_t *window, int cl);
void psg(int cl);
void write_windows(geometry_t *geom, unsigned long ***flag);
void write_window_map(geometry_t **windows, parameters_t *params);
void prints(double *A, char *fname, int k, int xlim, int ylim,
            double scale);
void get_local_ghost(geometry_t *geom, geometry_t *window, int n);
void windows_clear(hd_data_t *hd_data);
void reset_windows(hd_data_t *hd_data);
void reorder_cells(geometry_t *window, int *ncells, int *cells,
       int vc, int vcs, int *bot, int mode);
int get_local_sur(geometry_t *window, int cl, int ks);
void process_cell_mask(geometry_t *window, int *cells, int ncells);
void read_window_info(parameters_t *params, FILE *fp);
void write_site(geometry_t *window, double x, double y, char *tag);

/*------------------------------------------------------------------*/
/* Master routines                                                  */
/*------------------------------------------------------------------*/
master_t *master_alloc(void);
master_t *master_build(parameters_t *params, geometry_t *geom);
void master_free(master_t *master); /*UR-202 changed purpose and adjusted context */
void master_end(master_t *master);
void master_free_nwin(master_t *master);/*UR-202 changed  to reflect context */
void master_setghosts(geometry_t *geom, master_t *master, double **cellx, double **celly,
		      double **gridx, double **gridy, int *s2i, int *s2j);
void set_bdry_flags(geometry_t *geom, dump_data_t *dumpdata);
void compute_constants(parameters_t *params, geometry_t *geom,
                       master_t *master);
void pre_run_setup(master_t *master, geometry_t **window,
                   window_t **windat, win_priv_t **wincon);


/*------------------------------------------------------------------*/
/* Unstructured                                                     */
/*------------------------------------------------------------------*/
int zm1e(geometry_t *sgrid, int e);
int zp1e(geometry_t *sgrid, int e);
int e2e(geometry_t *sgrid, int e, int j);
void ete(geometry_t *sgrid, int e, int *e1, int *e2);
int jp(int j, int npe);
int jm(int j, int npe);
int jo(int j, int npe);
int jocw(geometry_t *geom, int c, int j);
void print_e(geometry_t *sgrid, mesh_t *m, int ed);
int read_grid_us(parameters_t *params, int cdfid);
int read_mesh_us(parameters_t *params);
void write_mesh_us(parameters_t *params, double *bathy, FILE *op, int mode);
void write_grid_us(parameters_t *params, geometry_t *geom);
void meshstruct_s(parameters_t *params, geometry_t *geom);
void meshstruct_us(parameters_t *params);
int get_mesh_obc(parameters_t *params, int **neic);
int get_mesh_obc_limit(parameters_t *params, int **neic);
void convert_mesh_input(parameters_t *params, mesh_t *mesh, double **xc, double **yc, 
			double *bathy, int wf);
void set_params_mesh(parameters_t *params, mesh_t *mesh, double **x, double **y);
mesh_t *mesh_init(parameters_t *params, int ns2, int npe);
void free_mesh(mesh_t *mesh);
int perimeter_mask(parameters_t *params, int **neic);
void write_mesh_desc(parameters_t *params, coamsh_t *cm, FILE *fp, int mode);
double edge_mean(geometry_t *window, double *a, int c);
double vertex_mean(geometry_t *window, double *a, int c);
void edge_centre(geometry_t *window, double *a, double *b, int mode);
void vertex_centre(geometry_t *window, double *a, double *b, int mode);
void mesh_reduce(parameters_t *params, double *bathy, double **xc, double **yc);
void interp_us(geometry_t *geom, double *vals, int nvec, int *vec, 
	       double *locx, double *locy, double *ret);
void interp_edge_us(geometry_t *geom, double *vals, int nvec, int *vec, 
		    double *locx, double *locy, double *ret, char *i_rule);
int is_index(geometry_t *geom, int cc, int c);

/*------------------------------------------------------------------*/
/* Open boundary routines                                           */
/*------------------------------------------------------------------*/
open_bdrys_t *OBC_alloc(void);
void OBC_build(open_bdrys_t **open, geometry_t *geom, geometry_t **window, int nwindows);
void get_OBC_conds(parameters_t *params, open_bdrys_t *open, FILE * fp, 
                   int n, tracer_info_t *tracers);
void init_OBC_conds(parameters_t *params, open_bdrys_t *open);
void get_obc_list(open_bdrys_t *open, FILE *fp, int n, char *key);
void convert_obc_list(parameters_t *params, open_bdrys_t *open, int n, geometry_t *geom, int *cmap);
void get_OBC_relax(parameters_t *params, open_bdrys_t *open, FILE *fp, int n);
void get_bdry_params(parameters_t *params, FILE *fp);
void set_OBC_cells(geometry_t *geom, open_bdrys_t *open, open_bdrys_t *io, 
		     int n, int nce1,
		     int nce2, int nz, unsigned long ***flag);
void bdry_custom_m(parameters_t *params, geometry_t *geom,
                   master_t *master);
void bdry_custom_w(geometry_t *geom, geometry_t **window);
void bdry_init_m(master_t *master);
void bdry_init_w(geometry_t *window, open_bdrys_t *open,
                 open_bdrys_t *gopen);
void bdry_custom_free(geometry_t *geom, master_t *master, geometry_t **window,
		      open_bdrys_t *open);
void bdry3d_m(geometry_t *geom, master_t *master);
void bdry2d_m(geometry_t *geom, master_t *master);
void copy_OBC_conds(open_bdrys_t *open, open_bdrys_t *io, int n,
                    tracer_info_t *tracer);
void bdry_data_copy(bdry_details_t *data, bdry_details_t *din);
void obc_loc(geometry_t *geom, int c, int *map, int *obc, int *oi1, int *oi2,
             int *nb2, int *nb3);
void obc_loc_tan(geometry_t *geom, int c, int *map, int *obc, int *oi1, int *oi2, 
		   int *nb2, int *nb3, unsigned long ***flag, unsigned long mask,
		   int i, int j, int k);
void bdry_tracer(geometry_t *window, window_t *windat, win_priv_t *wincon);
void set_OBC(geometry_t *window, window_t *windat, win_priv_t *wincon,
             open_bdrys_t *open, int sb, int eb, int *obc, int *oi1,
             int *oi2, int *cyc, double *vel, double *vel_t, double *vel_b,
	     int bcond, double dt, bdry_details_t *data,
	     double *transfer, int rlxn, int code);
void set_OBC_tr(int tn, geometry_t *window, window_t *windat,
                win_priv_t *wincon, open_bdrys_t *open);
void OBC_bgz_nograd(geometry_t *window);
void OBC_bgz_nogradb(geometry_t *window, open_bdrys_t *open, double *tr);

void bdry_eta(geometry_t *window, window_t *windat, win_priv_t *wincon);
void bdry_eval_u1_m(geometry_t *geom, master_t *master);
void bdry_eval_u2_m(geometry_t *geom, master_t *master);
void bdry_eval_tr_m(geometry_t *geom, master_t *master);
void bdry_eval_eta_m(geometry_t *geom, master_t *master);
void bdry_eval_u1av_m(geometry_t *geom, master_t *master);
void bdry_eval_u2av_m(geometry_t *geom, master_t *master);
void bdry_transfer_u1o(master_t *master, geometry_t **window,
                       window_t **windat);
void bdry_transfer_u1(master_t *master, geometry_t *window,
                      window_t *windat);
void bdry_transfer_u2o(master_t *master, geometry_t **window,
                       window_t **windat);
void bdry_transfer_u2(master_t *master, geometry_t *window,
                      window_t *windat);
void bdry_transfer_eta1(master_t *master, geometry_t **window,
                        window_t **windat);
void bdry_transfer_eta(master_t *master, geometry_t *window,
                       window_t *windat);
void bdry_transfer_u1av(master_t *master, geometry_t *window,
      window_t *windat);
void bdry_transfer_u2av(master_t *master, geometry_t *window,
      window_t *windat);
void bdry_transfer_tro(master_t *master, geometry_t **window,
                       window_t **windat);
void bdry_transfer_tr(master_t *master, geometry_t *window,
                      window_t *windat);
double bdry_value_w(geometry_t *window, window_t *windat,
                    win_priv_t *wincon, open_bdrys_t *open,
                    bdry_details_t *data, int c, int cc, double t);
double bdryval_t(geometry_t *window, window_t *windat, win_priv_t *wincon,
                 open_bdrys_t *open, bdry_details_t *data, int c, int c2,
                 double t);
void bcname(int code, char *name);
int sforce_relax(int bcond);
int cyclic_m1(geometry_t *geom, int codex, int codey, int *map, int c);
int cyclic_m2(geometry_t *geom, int codex, int codey, int *map, int c);
void set_sponge(geometry_t *window, double *u1vh, double dt, double *mask);
void set_sponge_c(geometry_t *window, double *AH, double dt);
void set_sponge_e(geometry_t *window, double *AH, double dt);
void reset_sponge_zone(geometry_t *window);
void set_sponge_cells(geometry_t *window);
void csr_tide_init(master_t *master, geometry_t **window);
void csr_tide_grid_init(master_t *master, geometry_t **window);
void csr_tide_grid_init_uv(master_t *master, geometry_t **window, int mode);
double csr_tide_eval(tidal_consts_t *tc, int cc, double jd);
void equ_tide_eval(geometry_t *window, window_t *windat, win_priv_t *wincon, double *equitide);
void custom_tide_init(master_t *master, geometry_t **window, int mode);
void custom_tide_grid_init(master_t *master, geometry_t **window);
void custom_tide_grid_init_uv(master_t *master, geometry_t **window, int mode);
double bc_tidal(geometry_t *window, window_t *windat, win_priv_t *wincon,
                tide_details_t *tide, int c);
void u1av_local(geometry_t *window, window_t *windat, win_priv_t *wincon,
		open_bdrys_t *open, double *val, int sb, int eb, int dir);
void eta_local(geometry_t *window, window_t *windat, win_priv_t *wincon,
	       open_bdrys_t *open, double *val, int sb, int eb);
void reset_sigma_OBC(geometry_t *window, window_t *windat, win_priv_t *wincon,
		     double *vel, 
		     double (*Hn) (geometry_t *, window_t *, win_priv_t *, int), 
		     int sb, int eb, int *obc);
void std_bdry(open_bdrys_t *open, int nfiles, char **file, 
	      tracer_info_t *tracers, double dt, double dt2d, char *bdrypath, int mode);
void bdry_reconfigure(master_t *master, geometry_t **window);
int bdry_reinit(geometry_t *geom, master_t *master, geometry_t **window, 
		open_bdrys_t *open, char *bstd);
void write_bdry(parameters_t *params, geometry_t *geom);
void std_bdry_init(open_bdrys_t *open, char *stdbdry, double dt, double dt2d,
		   char *bdrypath, tracer_info_t *tracer);
void upstrm(geometry_t *window,  window_t *windat, win_priv_t *wincon, open_bdrys_t *open,
	    double *newval, double *tr, int mode);
void average_OBC_corner(geometry_t *window, open_bdrys_t *open, double *var, int mode);
void reset_bdry_eta(geometry_t *window, window_t *windat, win_priv_t *wincon,
		    open_bdrys_t *open, double *eta);
void reset_obc_adjust(geometry_t *geom, double dt);

/*------------------------------------------------------------------*/
/* Distributed tranfer routines                                     */
/*------------------------------------------------------------------*/
void master_fill(master_t *master, geometry_t **window, window_t **windat,
                 win_priv_t **wincon);
void master_fill_ts(master_t *master, geometry_t **window, window_t **windat,
		    win_priv_t **wincon);
void master_fill_glider(master_t *master, geometry_t **window, window_t **windat,
			win_priv_t **wincon, ts_point_t *ts, double t);
void windat_fill(master_t *master, geometry_t *window, window_t *windat,
                 int nwindows, int mode);
void build_transfer_maps(geometry_t *geom, geometry_t **window, int wn,
			 int nwindows);
void win_data_fill_3d(master_t *master, geometry_t *window,
                      window_t *windat, int nwindows);
void win_data_refill_3d(master_t *master, geometry_t *window,
                        window_t *windat, int nwindows, int mode);
void win_data_empty_3d(master_t *master, geometry_t *window,
                       window_t *windat, int mode);
void win_data_fill_2d(master_t *master, geometry_t *window,
                      window_t *windat, int nwindows);
void win_data_refill_2d(master_t *master, geometry_t *window,
                        window_t *windat, int nwindows, int mode);
void win_data_empty_2d(master_t *master, geometry_t *window,
                       window_t *windat, int mode);
void window_reset(master_t *master, geometry_t *window, window_t *windat, 
		  win_priv_t *wincon, int mode);
void get_sloc(geometry_t *geom, int c);
void s2c_2d(geometry_t *geom, double *as, double **ac, int nx, int ny);
void v2c_2d(geometry_t *geom, double *as, double **ac, int nx, int ny);
void e2c_2d(geometry_t *geom, double *as, double **ac, int nx, int ny, int mode);
void s2c_3d(geometry_t *geom, double *as, double ***ac, int nx, int ny,
            int nz);
void e2c_3d(geometry_t *geom, double *as, double ***ac, int nx, int ny,
	     int nz, int mode);
void s2c_3d_e1(geometry_t *geom, double *as, double ***ac, int nx, int ny,
	       int nz);
void s2c_3d_e2(geometry_t *geom, double *as, double ***ac, int nx, int ny,
	       int nz);
void s2c_sed(geometry_t *geom, double **as, double ***ac, int nx, int ny,
       int nz);
void c2s_2d(geometry_t *geom, double *as, double **ac, int nx, int ny);
void c2s_3d(geometry_t *geom, double *as, double ***ac, int nx, int ny,
            int nz);
void c2e_2d(geometry_t *geom, double *as, double **ac, int nx, int ny, int mode);
void c2e_3d(geometry_t *geom, double *as, double ***ac, int nx, int ny,
            int nz, int mode);
void c2v_2d(geometry_t *geom,  double *as, double **ac,
            int nx, int ny);
void c2v_3d(geometry_t *geom,  double *as, double ***ac,
            int nx, int ny, int nz);
void s2m_3d(master_t *master, geometry_t *window, window_t *windat,
            win_priv_t *wincon);
void s2m_2d(master_t *master, geometry_t *window, window_t *windat);
void s2m_vel(double *Ag, double *Al, int *vec, int *evec, int nvec);
void s2m_flux(geometry_t *geom, geometry_t *window, double *Ag, double *Al,
              int *vec, int *evec, int nvec, int mode);
void set_precond(geometry_t *window, double *Al, int *vec, int nvec);
void pack_sparse(int *map, int mapsize, double *var, double *pack);
void unpack_sparse(int *map, int mapsize, double *var, double *pack, int oset);
void unpack_sparse3(int *map, int mapsize, double *var, double *pack, int oset);
void check_transfers(geometry_t *geom, geometry_t **window,
		     window_t **windat, win_priv_t **wincon, 
		     int nw, int mode);
void check_s2m();
int mpi_check_multi_windows_velocity(geometry_t *geom,master_t *master, geometry_t **window,
				     window_t **windat, int flag);
int mpi_check_multi_windows_Vz(geometry_t *geom,master_t *master, geometry_t **window,
			       window_t **windat);
int mpi_check_multi_windows_tracer(geometry_t *geom,master_t *master, geometry_t **window,
				     window_t **windat);
int mpi_check_multi_windows_sparse_arrays(geometry_t *geom);

 /*------------------------------------------------------------------*/
/* Transport routines                                               */
/*------------------------------------------------------------------*/
void trans_reset_init(parameters_t *params, master_t *master);
void trans_reset_end(master_t *master);
void transport_step(master_t *master, geometry_t **window,
		    window_t **windat, win_priv_t **wincon);
void transport_post(master_t *master, geometry_t **window,
		    window_t **windat, win_priv_t **wincon);
void read_tmode_params(FILE *fp, parameters_t *params);
void advect_diffuse_lag(geometry_t *window, window_t *windat,
			win_priv_t *wincon);
void calc_tmass(geometry_t *window, window_t *windat, win_priv_t *wincon,
		double *tmass, double *tsmass);
void get_weights(geometry_t *window, window_t *windat, win_priv_t *wincon);
void vel_center(geometry_t *window, window_t *windat, win_priv_t *wincon,
                double *nu, double *nv, double *nw);
void vel_center_w(geometry_t *window, window_t *windat,
		  win_priv_t *wincon, double *nw);
void semi_lagrange_c(geometry_t *window, window_t *windat,
                     win_priv_t *wincon);
void semi_lagrange_t(geometry_t *window, win_priv_t *wincon, geometry_t *tpg,
                     window_t *tpd, win_priv_t *tpc);
void semi_lagrange_tu(geometry_t *window, win_priv_t *wincon, geometry_t *tpg,
		      window_t *tpd, win_priv_t *tpc);
void semi_lagrange(geometry_t *window, window_t *windat,
                   win_priv_t *wincon, double *tr);
void streamline_atc(geometry_t *window, window_t *windat, win_priv_t *wincon, int c);
void semi_lagrange_atc(geometry_t *window, window_t *windat, win_priv_t *wincon, 
		       double *tr, double *ntr, int c);
double semi_lagrange_rtc(geometry_t *window, window_t *windat, win_priv_t *wincon, 
			 double *tr, int c);
void prep_semi_lagrange_atc(geometry_t *window);
void trans_data_check(master_t *master, geometry_t **window,
		      window_t **windat, win_priv_t **wincon);
void set_dzz(geometry_t *window, double *dzz);
int r2tij(geometry_t *window, int c, double p, double q, double *x, double *y);
void set_lmap(geometry_t *window, win_priv_t *wincon);
double int_val(geometry_t *window, int c, double *wgt, double *in);
double int_val_bl(geometry_t *window, int c, double *wgt, double *in);
double int_valo(geometry_t *window, double *wgt, int *lmap, double *in, int n);
void sl_check(geometry_t *window, int *c, double *cx, double *cy, double *cz);
void reset_Aij(geometry_t *window, window_t *windat, win_priv_t *wincon);
void get_verr(geometry_t *window, window_t *windat, win_priv_t *wincon);
void set_tracer_eta(geometry_t *window, window_t *windat, win_priv_t *wincon, 
		    double dt, double *oeta, double *neta);
void set_unit_eta(geometry_t *window, window_t *windat, win_priv_t *wincon, 
		  double dt, double *oeta, double *neta, double *Fx, double *Fy);
void tran_grid_init(geometry_t *window, window_t *windat, win_priv_t *wincon);
double vertex_weighted_tr(geometry_t *window, window_t *windat, 
			  win_priv_t *wincon, double *tr, int v);
int find_cell(geometry_t *window, int cd, double xs, double ys,
	      double *xi, double *yi);

/*------------------------------------------------------------------*/
/* Global geometry preprocessing routines                           */
/*------------------------------------------------------------------*/
hd_data_t *hd_init(FILE * prmfd);
int find_nearest(geometry_t *geom, double x, double y);
void alloc_geom(geometry_t *geom, unsigned long mask);
void alloc_geom_us(geometry_t *geom, unsigned long mask);
void geom_free(master_t *master, int mode);
void geom_free_us(master_t *master, geometry_t *geom, int mode);
void build_sparse_map(parameters_t *params);
void build_sparse_grid(parameters_t *params, geometry_t *sgrid, int nce1, 
		       int nce2, int nz, double **bathy, double *layers, 
		       int nobc, int *npts, int **iloc, int **jloc, int *type);
void build_sparse_grid_us(parameters_t *params, geometry_t *sgrid,
			  int nz, double *bathy, double *layers);
void alloc_sgrid(geometry_t *sgrid);
void clear_sgrid(geometry_t *sgrid);
void check_sparse(geometry_t *geom, int nwc, int ngc,
                  int nsc, int num_se, int num_ic, int num_dc,
                  int num_mc, unsigned long ***flg);
unsigned long se_oc(geometry_t *geom, int i, int j, int k,
                    int *xp1, int *xm1, int *yp1, int *ym1, int mode);
unsigned long di_ic(geometry_t *geom, int i, int j, int k,
                    int *xp1, int *xm1, int *yp1, int *ym1);
void point_geom(geometry_t *window, geometry_t *geom);
void point_geom_us(geometry_t *window, geometry_t *geom);
void make_flags(parameters_t *params, unsigned long ***flag,
                double **bathy, double *layers, int nce1, int nce2,
                int nz);
void u1_flags(unsigned long ***flag, int nce1, int nce2, int nz);
void u2_flags(unsigned long ***flag, int nce1, int nce2, int nz);
void corner_flags(unsigned long ***flag, int nce1, int nce2, int nz);
unsigned long b_flags(int nce1, int nce2, int nz, unsigned long ***flag,
                      int i, int j, int k, int *fi, int *fj,
                      int mode, int btype);
int set_kbot(int i, int j, unsigned long **flag, short **kbot,
             int btype, int nce1, int nce2);
void sigma_flags(int nfe1, int nfe2, int nz, unsigned long ***flg);
void make_cycled(geometry_t *geom, open_bdrys_t *open,
		 unsigned long ***flag);		 
int isvalidc(int c, int sz, char *warn);
void set_map_bdry(geometry_t *geom);
void get_filter(geometry_t *geom);

/*------------------------------------------------------------------*/
/* Diagnostic routines                                              */
/*------------------------------------------------------------------*/
void order (double *a, double *b);
void monitor(master_t *master, geometry_t **window, int mode);
void mass_diag(geometry_t *window, window_t *windat, win_priv_t *wincon);
void get_tend(geometry_t *window, int *vec, int eb, double *vel,
              double *ovel, double *tendency);
void get_tendv(geometry_t *window, int *vec, int eb,
	       double *tend, double *tend1, double *tend2);
void alerts(master_t *master);
void alerts_w(geometry_t *window, int mode);
void alerts_init(master_t *master, geometry_t **window);
void master_alert_fill(master_t *master, geometry_t *window,
           window_t *windat);
void write_run_setup(hd_data_t *hd_data);
void calc_cfl(geometry_t *window, window_t *windat, win_priv_t *wincon);
double int_wave_spd(geometry_t *window, window_t *windat,
                    win_priv_t *wincon, int c, int cc);
double int_wave_spd_cont_w(geometry_t *window, window_t *windat,
                           win_priv_t *wincon, int c, int cc);
double int_wave_spd_cont_m(master_t *master, int c, int cc);
double get_global_rms(double *a1, double *a2, int vc, int *vec);
double get_profile_rms(double *a1, double *a2, int cs, int *zm1);
void total_mass(geometry_t *window, window_t *windat, win_priv_t *wincon);
void steric(geometry_t *window, window_t *windat, win_priv_t *wincon);
void vorticity(geometry_t *window, window_t *windat, win_priv_t *wincon);
void calc_monotonic(geometry_t *window, window_t *windat, win_priv_t *wincon);
void diag_numbers(geometry_t *window, window_t *windat, win_priv_t *wincon);
void nor_vert_prof(geometry_t *window, window_t *windat, win_priv_t *wincon);
void ekman_pump_e1(geometry_t *window, window_t *windat, win_priv_t *wincon,
		   double *taus, double *taub);
void init_flushing(master_t *master, geometry_t **window,
                   win_priv_t **wincon);
void init_age(master_t *master, geometry_t **window, window_t **windat, 
                   win_priv_t **wincon);
void init_trperc(master_t *master, geometry_t **window, window_t **windat, 
                   win_priv_t **wincon);
void init_totals(master_t *master, geometry_t **window,
     win_priv_t **wincon);
void init_trans(master_t *master, geometry_t **window,
		win_priv_t **wincon);
void calc_flushing(master_t *master, geometry_t **window, window_t **windat, 
		   win_priv_t **wincon);
void reset_means(geometry_t *window, window_t *windat, win_priv_t *wincon,
                 int mode);
void reset_means_m_o(master_t *master);
void reset_means_m(master_t *master);
void init_means(master_t *master, parameters_t *params);
void get_means_w(geometry_t *window, window_t *windat, win_priv_t *wincon);
void get_means_m(geometry_t *geom, master_t *master,
                 dump_data_t *dumpdata);
void mld(geometry_t *window, window_t *windat, win_priv_t *wincon,
         int cs, int cb, double *tm, double *bm, int *css, int *ktb);
void mldk(geometry_t *window, window_t *windat, win_priv_t *wincon,
          int cs, int cb, double *tm, double *bm, int *css, int *ktb);
void check_nan(geometry_t *window, double *A, int as, int ae, int *vec,
               char *tag);
void check_max(geometry_t *window, double *A, int as, int ae, int *vec,
               char *tag);
double print_max(geometry_t *window, double *A, int as, int ae, int *vec,
     char *tag, FILE *fp, int *cl);
double get_max(geometry_t *window, double *A, int as, int ae,
         int *vec, int *cl);
double get_maxs(geometry_t *window, double *A, int as, int ae,
    int *vec, int *cl, double vmax);
void check_min(geometry_t *window, double *A, int as, int ae, int *vec,
               char *tag);
double get_min(geometry_t *window, double *A, int as, int ae,
         int *vec, int *cl);
double check_cont(geometry_t *window, int c);
int check_unstable(geometry_t *window, window_t *windat, win_priv_t *wincon,
		   int mode);
void debug_c(geometry_t *window, int var, int mode);
double percs(double *aa, int n, double pc);
void sorts(double *a, int n);
void sorts2(double *a, double *b, int *c, int n);
void orders(double *p, double *q);
void orders2(double *p, double *q, double *r, double *s, int *i, int *j);
int edge_sort_double (void const *x1, void const *x2);
int edge_sort_compare (void const *x1, void const *x2);
void perc_diag(geometry_t *window, window_t *windat, win_priv_t *wincon);
void calc_perc(FILE *fp);
void check_flux(geometry_t *window, window_t *windat, win_priv_t *wincon,
                int dw, int dc);
double is_night(double lat, double t);
void init_regions_g(master_t *master, parameters_t *params);
void init_regions_w(master_t *master, geometry_t **window);
void region_clean(geometry_t *win);
void region_transfer(master_t *master, geometry_t *window);
void region_mass(geometry_t *window, window_t *windat, win_priv_t *wincon);
void region_flux_coup(geometry_t *window, window_t *windat, win_priv_t *wincon,
		      double *Fx, double *Fz, double dt, int trn);
void region_flux_trans(geometry_t *window, window_t *windat, win_priv_t *wincon, int trn,
		       double *dtracer);
void region_mass_tr(geometry_t *window, window_t *windat, win_priv_t *wincon, 
		    double sign, int trn, int mode);
void region_volume_flux_coup(geometry_t *window, window_t *windat, win_priv_t *wincon,
			     double dt, double trem);
void region_volume_flux_trans(geometry_t *window, window_t *windat, win_priv_t *wincon);
void dump_regions(master_t *master);
void region_write(master_t *master, region_t *region);
void region_print(master_t *master, region_t *region);
void read_region_info(parameters_t *params, FILE *fp);
int read_regioni(master_t *master, char *rname, double *regionid);
int read_region_us(geometry_t *geom, char *name, double *regionid);
double con_velmax(double vel, double velmax);
double con_mean(geometry_t *window, double *vel, double velmax, int c, int mode);
double con_median(geometry_t *window, double *vel, double velmax, int c, int mode);
double con_med(geometry_t *window, double *vel, int c);
double decorr(geometry_t *window, double *a, int sz, int c, int mode, double scale);
void calc_decorr(geometry_t *window, double *a, double *dex, int sz, int mode, 
		 double scale);
void calc_dhd(geometry_t *window, window_t *windat, win_priv_t *wincon);
double buoyancy_frequency2_m(master_t * master, geometry_t *geom, double *dens, int c);

/*------------------------------------------------------------------*/
/* Forcing routines                                                 */
/*------------------------------------------------------------------*/
int yrday(int year, int mon, int day);
void forcings_init(master_t *master);
void forcings_end();
void frc_ts_eval_grid(master_t *master, double t, timeseries_t *ts, int id,
                      double *p, double conv);
void frc_ts_eval_grid_mult(master_t *master, double t, timeseries_t **ts, int *id,
			   int ntsfiles, double *p, double conv);
timeseries_t *frc_read_cell_ts_o(master_t *master, char *key,
         char *varname, char *varunit, double *dt,
         int *id, double **p);
timeseries_t *frc_read_cell_ts(master_t *master, char *fname, double i_dt,
             char *varname, char *varunit, double *dt,
             int *id, double **p);
timeseries_t *frcw_read_cell_ts(master_t *master, char *fname, double i_dt,
             char *varname, char *varunit, double *dt,
             int *id, double **p);
timeseries_t **frc_read_cell_ts_mult(master_t *master, char *fname, double i_dt,
				     char *varname, char *varunit, double *dt,
				     int *id, int *ntsfiles, double **p, int quitmode);
double frc_get_input_dt(double dt, double wdt, char *name);
/* Heatflux                                                         */
void calc_heatf(geometry_t *window, window_t *windat, win_priv_t *wincon);
void calc_saltf(geometry_t *window, window_t *windat, win_priv_t *wincon);
double large_pond(double sst, double ta, double ws, double zref,
		  double *c_e, double *c_h);
double kitiagorodskii(double ws, double zref, double *c_e, double *c_h);
double bunker(double sst, double ta, double ws, double zref,
	      double *c_e, double *c_h);
double kondo(double sst, double ta, double wspd, double qs, double qq,
	     double zref, double *c_e, double *c_h);
/* Wind                                                             */
void windstress(double *wx, double *wy, double v0, double v1,
                double c0, double c1);
void stresswind(double *wx, double *wy, double v0, double v1,
                double c0, double c1);
int wind_init(sched_event_t *event);
double wind_event(sched_event_t *event, double t);
void wind_cleanup(sched_event_t *event, double t);

int patm_init(sched_event_t *event);
double patm_event(sched_event_t *event, double t);
void patm_cleanup(sched_event_t *event, double t);
int patm_init_single(sched_event_t *event);
double patm_event_single(sched_event_t *event, double t);
void patm_cleanup_single(sched_event_t *event, double t);

int precip_init(sched_event_t *event);
double precip_event(sched_event_t *event, double t);
void precip_cleanup(sched_event_t *event, double t);

int evap_init(sched_event_t *event);
double evap_event(sched_event_t *event, double t);
void evap_cleanup(sched_event_t *event, double t);

int airtemp_init(sched_event_t *event);
double airtemp_event(sched_event_t *event, double t);
void airtemp_cleanup(sched_event_t *event, double t);

int rh_init(sched_event_t *event);
double rh_event(sched_event_t *event, double t);
void rh_cleanup(sched_event_t *event, double t);

int wetbulb_init(sched_event_t *event);
double wetbulb_event(sched_event_t *event, double t);
void wetbulb_cleanup(sched_event_t *event, double t);

int cloud_init(sched_event_t *event);
double cloud_event(sched_event_t *event, double t);
void cloud_cleanup(sched_event_t *event, double t);

int swr_init(sched_event_t *event);
double swr_event(sched_event_t *event, double t);
void swr_cleanup(sched_event_t *event, double t);

int light_init(sched_event_t *event);
double light_event(sched_event_t *event, double t);
void light_cleanup(sched_event_t *event, double t);

int webf_init(sched_event_t *event);
double webf_event(sched_event_t *event, double t);
void webf_cleanup(sched_event_t *event, double t);
void webf_write_setup(FILE *fp, sched_event_t *event);

int heatf_init(sched_event_t *event);
double heatf_event(sched_event_t *event, double t);
void heatf_cleanup(sched_event_t *event, double t);

int longout_init(sched_event_t *event);
double longout_event(sched_event_t *event, double t);
void longout_cleanup(sched_event_t *event, double t);

int longin_init(sched_event_t *event);
double longin_event(sched_event_t *event, double t);
void longin_cleanup(sched_event_t *event, double t);

int sheat_init(sched_event_t *event);
double sheat_event(sched_event_t *event, double t);
void sheat_cleanup(sched_event_t *event, double t);

int lheat_init(sched_event_t *event);
double lheat_event(sched_event_t *event, double t);
void lheat_cleanup(sched_event_t *event, double t);

int swr_mean_init(sched_event_t *event);
double swr_mean_event(sched_event_t *event, double t);
void swr_mean_cleanup(sched_event_t *event, double t);

int hftemp_init(sched_event_t *event);
double hftemp_event(sched_event_t *event, double t);
void hftemp_cleanup(sched_event_t *event, double t);

int storm_init(sched_event_t *event);
double storm_event(sched_event_t *event, double time);
void storm_cleanup(sched_event_t *event, double t);

void crash_recovery_init(hd_data_t *hd_data);
double crash_event(sched_event_t *event, double t);

int regulate_init(sched_event_t *event);
double regulate_event(sched_event_t *event, double time);
void regulate_cleanup(sched_event_t *event, double t);

void swr_params_init(master_t *master, geometry_t **window);
double swr_params_event(geometry_t *window, window_t *windat, win_priv_t *wincon, int n);

/*------------------------------------------------------------------*/
/* Data Assimilation routines                                       */
/*------------------------------------------------------------------*/
#ifdef HAVE_DA
void dassim_init(master_t *master);
void dassim_post_init(master_t *master);
void dassim_end(void);
void dassim_run_setup(FILE *fp, sched_event_t *event);
#endif

/*------------------------------------------------------------------*/
/* Grid refinement subroutines                                      */
/*------------------------------------------------------------------*/
void build_zoom_maps(geometry_t *geom, geometry_t **window);
void set_zoom_m2d(geometry_t *window);
void hgrid_zoom(geometry_t *geom, geometry_t *window);
void hgrid_zoom_p(geometry_t *geom, geometry_t **window);
void hgrid_zoom_m(geometry_t *geom, geometry_t **window, dump_data_t *dumpdata);
void consts_zoom(master_t *master, geometry_t *geom,
		 geometry_t *window, win_priv_t *wincon);

void zflux_e1(geometry_t *geom, geometry_t *window, double *Ag, double *Al,
              int *vec, int *evec, int nvec);
void zflux_e2(geometry_t *geom, geometry_t *window, double *Ag, double *Al,
              int *vec, int *evec, int nvec);
void zvel_filter(geometry_t *geom, geometry_t *window, double *Ag, double *Al,
		 double zs, int *vec, int *evec, int nvec);
void zoom_shift(geometry_t *geom, geometry_t *window);
void set_intrp(geometry_t *geom, geometry_t *window);
void geom_interp(geometry_t *geom, geometry_t **window);
void geom_interp_e1(geometry_t *geom, geometry_t **window);
void geom_interp_e2(geometry_t *geom, geometry_t **window);
void geom_interp_e1(geometry_t *geom, geometry_t **window);
void geom_interp_e2(geometry_t *geom, geometry_t **window);
void global_interp(geometry_t *geom, double *A, int nz);
void global_interp_e1(geometry_t *geom, double *A, int nz);
void global_interp_e2(geometry_t *geom, double *A, int nz);
void global_sinterp_e1(geometry_t *geom, double *A, int nzp, int nzc);
void global_sinterp_e2(geometry_t *geom, double *A, int nzp, int nzc);
void global_interp_nn(geometry_t *geom, double *A, int nz);
void global_interp_e1_nn(geometry_t *geom, double *A, int nz);
void global_interp_e2_nn(geometry_t *geom, double *A, int nz);
void global_eta(geometry_t *geom, double dt2d);
void global_interp_flux(master_t *master, geometry_t *geom, int ic,
      double dt2d);
void set_zoom_flux(geometry_t *window, window_t *windat, win_priv_t *wincon);
void global_interp_depth(master_t *master, geometry_t *geom);
void global_interp_vel2d(master_t *master, geometry_t *geom);

void global_interp_pre2d(master_t *master, geometry_t *geom, int ic,
			  double dt2d);

void global_interp_post2d(master_t *master, geometry_t *geom);
void global_interp_pre(master_t *master, geometry_t *geom);
void global_interp_post(master_t *master, geometry_t *geom);
void interp_adjust_e1(master_t *master, geometry_t *geom);
void interp_adjust_e2(master_t *master, geometry_t *geom);
void interp_dz_at_u1(master_t *master, geometry_t *geom);
void interp_dz_at_u2(master_t *master, geometry_t *geom);
void glob_fill3D(geometry_t *geom, master_t *master);
void glob_fill2D(geometry_t *geom, master_t *master);
void fill_zoom_cell(geometry_t *geom, double *A, int c, int ze1, int ze2);
int zoom_search(geometry_t *geom, geometry_t *window, int gc,
                int lc, int *xmap, int *ymap, int *ii, int *jj);
void interp_zoom(geometry_t *window, double *A, int nz);
void interp_flux(geometry_t *window, double *A, int nz, int mode);
void interp_zoom1(geometry_t *window, double *A, int *vec, int nvec);
void interp_flux1(geometry_t *window, double *A, int *vec, int nvec,
                  int mode);
void set_window_codes(geometry_t *geom, geometry_t *window, int *vec,
                      int nvec, int zc);
void reset_zoom_ghosts(geometry_t *geom, geometry_t *window, int *vec,
                       int nvec);
void set_sponge_zoom(geometry_t *window, window_t *windat, win_priv_t *wincon,
         double *u1vh, double *u2vh, int nzbe1, int nzbe2);
void set_flux_scale(geometry_t *geom, geometry_t *window);
void set_blend_zones(master_t *master, geometry_t *window, win_priv_t *wincon);

/*------------------------------------------------------------------*/
/* Tracer subroutines                                               */
/*------------------------------------------------------------------*/
int tracer_write(parameters_t *params, FILE *op);
void tracer_write_3d(master_t *master, FILE *op, tracer_info_t *tracer, int n);
void tracer_write_2d(master_t *master, FILE *op, tracer_info_t *tracer, int n);
void tracer_step(master_t *master, geometry_t **window, window_t **windat,
                 win_priv_t **wincon, int nwindows);
void tracer_step_3d(geometry_t *window, window_t *windat, win_priv_t *wincon);
void tracer_step_2d(geometry_t *window, window_t *windat, win_priv_t *wincon);
void tracer_step_init(master_t *master);
void create_tracer_3d(parameters_t *params);
void init_tracer_2d(parameters_t *params, master_t *master);
void init_tracer_3d(parameters_t *params, master_t *master);
void init_tracer_sed(parameters_t *params, master_t *master);
void tracer_relax_init(master_t *master);
void tracer_relax_end(master_t *master);
void tr_relax_event_c(sched_event_t *event, double relax_rate,
          double dt, double t, int c);
void tracer_reset_init(master_t *master);
void tracer_reset_end(master_t *master);
void tracer_reset2d_init(master_t *master);
void tracer_reset2d_end(master_t *master);
void tracer_dhw_init(master_t *master);
void set_lateral_BC_tr(double **tr, int ntr, int sgbpt, int *bpt,
                       int *bin);
int advect_diffuse(geometry_t *window, window_t *windat, win_priv_t *wincon);
int advect_diffuse_2d(geometry_t *window, window_t *windat, win_priv_t *wincon);
void advect(geometry_t *window, window_t *windat, win_priv_t *wincon,
            double *tr, double *Fx, double *Fz, 
	    double *tr_mod, double *tr_mod_x, double *tr_mod_z, double dtu, int trasc);
void diffuse(geometry_t *window, window_t *windat, win_priv_t *wincon);
void hor_diffuse(geometry_t *window, window_t *windat, win_priv_t *wincon,
                 double *tr, double *Fx);
void hor_diffuse_2d(geometry_t *window, window_t *windat,
                    win_priv_t *wincon, double *tr, double *Fx);
void vert_diffuse_3d(geometry_t *window, window_t *windat,
                     win_priv_t *wincon);
void vert_diffuse_2d(geometry_t *window, window_t *windat,
                     win_priv_t *wincon);
void implicit_vdiff_tr(geometry_t *window, window_t *windat,
                       win_priv_t *wincon, double *var, double *Kz,
                       double *dzcell, double *dzface, double *fb,
                       double *ft, int *ctp, int *cbt, int vcs,
                       double *Splus, double *Sminus, double *scale,
                       double *C, double *Cp1, double *Cm1);
void implicit_vdiff_at_cc(geometry_t *window, window_t *windat,
			  win_priv_t *wincon, double *var, double *Kz,
			  double *dzcell, double *dzface, double *fb,
			  double *ft, int *ctp, int *cbt, int cc,
			  double *Splus, double *Sminus, double *scale,
			  double *C, double *Cp1, double *Cm1, double dt);
void mode2d_tracer_init(geometry_t *window, window_t *windat,
                        win_priv_t *wincon);
void tr_diag_reset_w(geometry_t *window, window_t *windat,
                     win_priv_t *wincon);
void tr_diag_reset_m(master_t *master);
void set_sbc(geometry_t *window, window_t *windat, win_priv_t *wincon);
void calc_swr(geometry_t *window, window_t *windat, win_priv_t *wincon,
	      double *So, int cci);
void surf_relax(geometry_t *window, window_t *windat, win_priv_t *wincon);
double surf_conc(geometry_t *window, window_t *windat, win_priv_t *wincon,
                 double *tr, double osubeta, double subeta, int c, int c2, int cc,
                 double fcbot, double *dtracer);
double surf_eta(geometry_t *window, window_t *windat, win_priv_t *wincon,
                 int c, int c2, int cc, double fcbot, double oldeta, double *dtracer);
int get_icell(geometry_t *geom, geometry_t *window, int cl, int cg, int ci);
void set_map_t(geometry_t *window);
void reset_map_t(geometry_t *window);
void reset_map_t_all(geometry_t *window);
void set_map_inside(geometry_t *window);
void set_map_l(geometry_t *window);
void ultimate_filter(double *cr, double *Fx, double *F, int e, 
		    int c, double trp, int cm1, double trm);
void set_dz(geometry_t *window, window_t *windat, win_priv_t *wincon);
void order1(geometry_t *window, window_t *windat, win_priv_t *wincon,
            double *tr, double *Fx, double *Fz, double dtu);
void order2(geometry_t *window, window_t *windat, win_priv_t *wincon,
            double *tr, double *Fx, double *Fz);
void order4(geometry_t *window, window_t *windat, win_priv_t *wincon,
            double *tr, double *Fx, double *Fz);
void order3us(geometry_t *window, window_t *windat, win_priv_t *wincon,
	      double *tr, double *Fx, double *Fz);
void order4us(geometry_t *window, window_t *windat, win_priv_t *wincon,
	      double *tr, double *Fx, double *Fz);
void order4_do(double *F, double *tr, int e, int cp1, int cp2, int cm1, int cm2);
void quickest(geometry_t *window, window_t *windat, win_priv_t *wincon,
              double *tr, double *Fx, double *Fz);
void quickest_uniform(geometry_t *window, window_t *windat,
                      win_priv_t *wincon, double *tr, double *Fx,
                      double *Fz);
void van_leer(geometry_t *window, window_t *windat, win_priv_t *wincon,
              double *tr, double *Fx, double *Fz);
void van_leer_do(double *F, double *tr, double *vel, double *cn, 
		 int e, int cp1, int cp2, int cm1, int cm2);
void van_leer_tr(double *F, double *tr, double *vel, double *cn, 
		 int e, double tp1, double tp2, double tm1, double tm2);
int intersect(geometry_t *window, int c, int j, double xn, double yn,
	      double xs, double ys, double *xi, double *yi);
void ffsl_init(geometry_t *window, window_t *windat, win_priv_t *wincon);
void ffsl_trans_prep(geometry_t *window, window_t *windat, win_priv_t *wincon);
void ffsl_trans_post(geometry_t *window, window_t *windat, win_priv_t *wincon);
void prep_ff_sl(geometry_t *window, window_t *windat, win_priv_t *wincon,
	       double dt);
void ff_sl_do_vert(double *vel, double dt, double *h, int ss,
	     int se, int *sdo, int *fmap, int *bmap, int *cl, double *crf);
void ff_sl_van_leer(double *F, double *tr, double *vel, double *crf,
		   int *cl, int ss, int se, int *sdo, int *fmap, int *bmap);
void ffsl_do(geometry_t *window, window_t *windat, win_priv_t *wincon,
	     double *tr, double *Fx, double *Fz, double dtu);
void ffsl_don(geometry_t *window, window_t *windat, win_priv_t *wincon,
	      double *tr, double *Fx, double *Fz, double dtu);
double hd_trans_interp(geometry_t *window, GRID_SPECS **gs, double x, double y, double z, 
		       int c, int co, int vid);
void clip_ffsl(geometry_t *window, window_t *windat, win_priv_t *wincon,
	       double *tr);
void order2_upwind(geometry_t *window, window_t *windat, win_priv_t *wincon,
		   double *tr, double *Fx, double *Fz);
void order2_upwind_do(double *F, double *tr, double *vel, double *cn, 
		      int e, int cp1, int cp2, int cm1, int cm2);
void set_multidt_t(geometry_t *window, window_t *windat, double *tr,
                   double trem, int tn);
void build_advect_weights(geometry_t *window, window_t *windat, win_priv_t *wincon);
void get_advect_metrics(geometry_t *window, window_t *windat, win_priv_t *wincon, 
			int e);
void build_quadratic_weights(geometry_t *window, window_t *windat, win_priv_t *wincon);
void get_quadratic_metrics(geometry_t *window, window_t *windat, win_priv_t *wincon, 
			   int e);
double get_quadratic_value(geometry_t *window, window_t *windat, win_priv_t *wincon, 
			   double *tr, int e, int j);
void build_linear_weights(geometry_t *window, window_t *windat, win_priv_t *wincon);
void get_linear_metrics(geometry_t *window, window_t *windat, win_priv_t *wincon, int co,
			double **xc, double **yc, double **zc);
void get_linear_limit(geometry_t *window, window_t *windat, win_priv_t *wincon, double *tr);
double get_linear_value(geometry_t *window, window_t *windat, win_priv_t *wincon, double *tr, int co,
			double x, double y, double z);
void get_local_bounds(geometry_t *window, window_t *windat, win_priv_t *wincon, double tr, double *trv,
		      int e, int co, double *smin, double *smax, int mode);
void get_wet_cells(geometry_t *window, win_priv_t *wincon);
void calc_courant(geometry_t *window, window_t *windat, win_priv_t *wincon,
                  double dt);
void get_surface(geometry_t *window, window_t *windat, win_priv_t *wincon,
                 double *tr, double *Fz, double *dtracer, double mintr,
                 double maxtr);
void update_flux(geometry_t *window, window_t *windat, win_priv_t *wincon,
                 double *Fx, double *Fz, double dtu);
void get_dtracer(geometry_t *window, window_t *windat, win_priv_t *wincon,
                 double *dtracer, double *Fx, double *Fy, double dtu);
void update_tracer(geometry_t *window, window_t *windat,
                   win_priv_t *wincon, double *tr, double *dtracer,
                   double *Fz, int ksf, int slf, double mintr,
                   double maxtr);
void set_flux_3d(geometry_t *window, window_t *windat, win_priv_t *wincon,
                 int mode);
void tr_bounds(geometry_t *window, window_t *windat, win_priv_t *wincon);
void tracer_decay(geometry_t *window, window_t *windat, win_priv_t *wincon,
                  double decay, double *tracer, tracer_info_t *tr);
void vel_center(geometry_t *window, window_t *windat, win_priv_t *wincon,
                double *nu, double *nv, double *nw);
void load_tracer_step_3d(parameters_t *params, master_t *master,
                         FILE * fp);
void load_tracer_step_2d(parameters_t *params, master_t *master,
                         FILE * fp);
void load_wc_tracer_name(master_t *master, FILE * fp, char *trname, int dim);
void temp_salt_init(parameters_t *params, master_t *master);
void set_multidt_tr_2d(geometry_t *window, window_t *windat, double trem,
                       double *vel, double *as, double *ae);
void set_lateral_BC_rad(geometry_t *window, window_t *windat, 
			win_priv_t *wincon);
int get_inc(char *buf);
void scale_tracer(tracer_info_t *trinfo, geometry_t *window, double **tr, int n);
void tracer_re_read(tracer_info_t *tr, FILE *fpt, int trinfo_type);
void calc_flux(geometry_t *window, window_t *windat, win_priv_t *wincon,
	       double *fluxe1, double *fluxw, double dt,
	       int dir1, int dir2);
void calc_flux_old(geometry_t *window, window_t *windat, win_priv_t *wincon, double dt);
void check_bounds(char *text, geometry_t *window, window_t *windat, 
		  win_priv_t *wincon);

/*------------------------------------------------------------------*/
/* 2D velocity subroutines                                          */
/*------------------------------------------------------------------*/
void mode2d_step(geometry_t *geom, master_t *master, geometry_t **window,
                 window_t **windat, win_priv_t **wincon, int nwindows);
void vel2D_lbc(double *vel, int sgbpt, int sgn, int *bpt,
               int *bin, double slip);
void set_lateral_BC_vel(double *vel, int sgbpt, int *bpt, int *bin);
void set_lateral_bc_eta(double *eta, int sgbpt, int *bpt, int *bin,
                        int *bin2, int mode);
void vel_u1av_update(geometry_t *window, window_t *windat,
                     win_priv_t *wincon);
void vel_u1av_update_seq(geometry_t *window, window_t *windat,
			 win_priv_t *wincon);
void set_map_e1(geometry_t *window);
void advect_u1_2d(geometry_t *window, window_t *windat, win_priv_t *wincon);
void advect_u1_2d_ang_adv_f(geometry_t *window, window_t *windat,
			    win_priv_t *wincon);
void advect_u1_2d_ang_adv_b(geometry_t *window, window_t *windat,
          win_priv_t *wincon);
void advect_u1_2d_ang_flux_f(geometry_t *window, window_t *windat,
           win_priv_t *wincon);
void advect_u1_2d_ang_flux_b(geometry_t *window, window_t *windat,
           win_priv_t *wincon);
void hvisc_setup_2d(geometry_t *window, window_t *windat,
		    win_priv_t *wincon, double *tzp);
void hvisc_setup_2d_old(geometry_t *window, window_t *windat,
			win_priv_t *wincon, double *tzp);
void hvisc_u1_2d(geometry_t *window, window_t *windat, win_priv_t *wincon,
                 double *tzp);
void hvisc_u1_2dus(geometry_t *window, window_t *windat, win_priv_t *wincon, 
		   double *tzp);
void hvisc_u1_2dusb(geometry_t *window, window_t *windat, win_priv_t *wincon, 
		    double *tzp);
void hvisc_u1_2d_simple(geometry_t *window, window_t *windat,
                        win_priv_t *wincon, double *tzp);
void hvisc_2d_null(geometry_t *window, window_t *windat, win_priv_t *wincon,
		   double *tzp);
void mdxs(geometry_t *window, win_priv_t *wincon);
void mdys(geometry_t *window, win_priv_t *wincon);
double mdxns(geometry_t *window, window_t *windat, win_priv_t *wincon,
             int c);
double mdxbs(geometry_t *window, window_t *windat, win_priv_t *wincon,
             int c);
double mdyns(geometry_t *window, window_t *windat, win_priv_t *wincon,
             int c);
double mdybs(geometry_t *window, window_t *windat, win_priv_t *wincon,
             int c);
void set_map_e2(geometry_t *window);
void bdry_u1_2d(geometry_t *window, window_t *windat, win_priv_t *wincon);
void advect_u2_2d(geometry_t *window, window_t *windat, win_priv_t *wincon);
void advect_u2_2d_ang_adv_f(geometry_t *window, window_t *windat,
          win_priv_t *wincon);
void advect_u2_2d_ang_adv_b(geometry_t *window, window_t *windat,
          win_priv_t *wincon);
void advect_u2_2d_ang_flux_f(geometry_t *window, window_t *windat,
           win_priv_t *wincon);
void advect_u2_2d_ang_flux_b(geometry_t *window, window_t *windat,
           win_priv_t *wincon);
void semi_lagrange_u1av(geometry_t *window, window_t *windat, win_priv_t *wincon);
void semi_lagrange_u2av(geometry_t *window, window_t *windat, win_priv_t *wincon);
void hvisc_u2_2d(geometry_t *window, window_t *windat, win_priv_t *wincon,
                 double *tzp);
void hvisc_u2_2d_simple(geometry_t *window, window_t *windat,
                        win_priv_t *wincon, double *tzp);
void bdry_u2_2d(geometry_t *window, window_t *windat, win_priv_t *wincon);
void store_old_bdry_values(double *F, open_bdrys_t *open,
                           bdry_old_values_t *old, int sb, int eb,
                           int *obc, int *oi1, int *oi2);
void eta_step(geometry_t *window, window_t *windat, win_priv_t *wincon,
              int ic);
void set_map_eta(geometry_t *window);
void asselin(geometry_t *window, window_t *windat, win_priv_t *wincon);
void leapfrog_update_2d(geometry_t *window, window_t *windat,
                        win_priv_t *wincon);
void extract_velocity_2d(geometry_t *window, window_t *windat,
                         win_priv_t *wincon);
void extract_u2(geometry_t *window, window_t *windat, win_priv_t *wincon);
void get_depths(geometry_t *window, window_t *windat, win_priv_t *wincon);
int eta_relax_init(sched_event_t *event);
double eta_relax_event(sched_event_t *event, double t);
void eta_relax_cleanup(sched_event_t *event, double t);
void do_eta_relax(geometry_t *window, int c, double relax_rate, int tidef);
void do_eta_increment(geometry_t *window);
int vel_relax_init(sched_event_t *event);
double vel_relax_event(sched_event_t *event, double t);
void vel_relax_cleanup(sched_event_t *event, double t);
void do_vel_relax(geometry_t *window, window_t *windat, win_priv_t *wincon, int mode);
void do_vel_increment_2d(geometry_t *window);
void vint_3d(geometry_t *window, window_t *windat, win_priv_t *wincon);

/*------------------------------------------------------------------*/
/* 3D velocity subroutines                                          */
/*------------------------------------------------------------------*/
void mode3d_step(geometry_t *geom, master_t *master, geometry_t **window,
                 window_t **windat, win_priv_t **wincon, int nwindows);
void mode3d_post(geometry_t *geom, master_t *master, geometry_t **window,
                 window_t **windat, win_priv_t **wincon, int nwindows);
void mode3d_prep(geometry_t *geom, master_t *master, geometry_t **window,
                 window_t **windat, win_priv_t **wincon, int nwindows);
void vel_u1_update(geometry_t *window, window_t *windat,
                   win_priv_t *wincon);
void vel_u2_update(geometry_t *window, window_t *windat,
                   win_priv_t *wincon);
void set_dz_at_u1(geometry_t *window, window_t *windat,
                  win_priv_t *wincon);
void set_thin_u1(geometry_t *window, window_t *windat, win_priv_t *wincon);
void cells2process_e1(geometry_t *window, window_t *windat,
                      win_priv_t *wincon);
void cells2process_e2(geometry_t *window, window_t *windat,
                      win_priv_t *wincon);
void linear_bdry_cell(geometry_t *window, window_t *windat, win_priv_t *wincon,
          int *ctp, int nctp, int *mask, int *bottom);
void set_dz_at_u2(geometry_t *window, window_t *windat,
                  win_priv_t *wincon);
void set_thin_u2(geometry_t *window, window_t *windat, win_priv_t *wincon);
void set_surf_cond(geometry_t *window, double *a, double val, int *vec, int nvec);
void set_surf_cells(geometry_t *window,  window_t *windat, 
		    win_priv_t *wincon, int mode);
int advect_u1_3d(geometry_t *window, window_t *windat, win_priv_t *wincon);                  
int nonlin_coriolis_3d(geometry_t *window, window_t *windat, win_priv_t *wincon);
int nonlin_coriolis_2d(geometry_t *window, window_t *windat, win_priv_t *wincon);
double pv_energy_neutral(window_t *windat, int e, int eoe, int v1, int v2);
double pv_enstrophy_conserve(window_t *windat, int e, int eoe, int v1, int v2);
double pv_enstrophy_dissipate(window_t *windat, int e, int eoe, int v1, int v2);
int advect_u2_3d(geometry_t *window, window_t *windat, win_priv_t *wincon);
void advect_u1_3d_ang_adv_f(geometry_t *window, window_t *windat,
			    win_priv_t *wincon);
void advect_u2_3d_ang_adv_f(geometry_t *window, window_t *windat,
			    win_priv_t *wincon);
void advect_u1_3d_ang_flux_f(geometry_t *window, window_t *windat,
			     win_priv_t *wincon);
void advect_u2_3d_ang_flux_f(geometry_t *window, window_t *windat,
			     win_priv_t *wincon);
void advect_u2_3d_ang_adv_b(geometry_t *window, window_t *windat,
          win_priv_t *wincon);
void advect_u1_3d_ang_adv_b(geometry_t *window, window_t *windat,
          win_priv_t *wincon);
void advect_u2_3d_ang_flux_b(geometry_t *window, window_t *windat,
           win_priv_t *wincon);
void advect_u1_3d_ang_flux_b(geometry_t *window, window_t *windat,
           win_priv_t *wincon);
void semi_lagrange_u1(geometry_t *window, window_t *windat, win_priv_t *wincon,
		      int *cells, int vcs, int vc);
void semi_lagrange_u2(geometry_t *window, window_t *windat, win_priv_t *wincon,
		      int *cells, int vcs, int vc);
void hvisc_init(master_t *master, win_priv_t **wincon);
void hvisc_setup(geometry_t *window, window_t *windat, win_priv_t *wincon);
void hvisc_setup_old(geometry_t *window, window_t *windat, win_priv_t *wincon);
void hvisc_setup_pre(geometry_t *window, window_t *windat, win_priv_t *wincon);
void hvisc_setup_pre2d(geometry_t *window, window_t *windat, win_priv_t *wincon);
void hvisc_u1_3d(geometry_t *window, window_t *windat, win_priv_t *wincon);
void hvisc_u1_3dus(geometry_t *window, window_t *windat, win_priv_t *wincon);
void hvisc_u1_3dusb(geometry_t *window, window_t *windat, win_priv_t *wincon);
void hvisc_u1_3d_simple(geometry_t *window, window_t *windat,
                        win_priv_t *wincon);
void hvisc_null(geometry_t *window, window_t *windat, win_priv_t *wincon);
void reset_hdiff(geometry_t *window, int cl);
void set_hdiff(geometry_t *window, double *AH, double AH0);
void reset_hor_diff(master_t *master, double u1vh, int flag);
int *stencil(geometry_t *window, int cl, int *size, int type, int edge);
void pressure_u1(geometry_t *window, window_t *windat, win_priv_t *wincon);
void coriolis_u1(geometry_t *window, window_t *windat, win_priv_t *wincon);
void stokes_u1(geometry_t *window, window_t *windat, win_priv_t *wincon);
int vdiff_u1(geometry_t *window, window_t *windat, win_priv_t *wincon);
void bdry_u1_3d(geometry_t *window, window_t *windat, win_priv_t *wincon);
void extract_u1_3d(geometry_t *window, window_t *windat,
                   win_priv_t *wincon);
void leapfrog_update_3d(geometry_t *window, window_t *windat,
                        win_priv_t *wincon);
void init_sigma(geometry_t *window, window_t *windat, win_priv_t *wincon);
void velocity_adjust(geometry_t *window, window_t *windat,
                     win_priv_t *wincon);
void set_new_cells_u1(geometry_t *window, window_t *windat,
                      win_priv_t *wincon);
void vel_w_update(geometry_t *window, window_t *windat,
                  win_priv_t *wincon);
void vel_w_trans(geometry_t *window, window_t *windat, win_priv_t *wincon);
void ff_sl_w_update(geometry_t *window, window_t *windat, win_priv_t *wincon);
void init_wvel_bounds(geometry_t *geom, master_t *master);
void vel_w_bounds(geometry_t *window, window_t *windat, win_priv_t *wincon);
void vel_w_bounds_tran(geometry_t *window, window_t *windat, win_priv_t *wincon);
void bdry_w(geometry_t *window, window_t *windat, win_priv_t *wincon);
int implicit_vdiff(geometry_t *window, window_t *windat,
		   win_priv_t *wincon, double *var, double *nvar,
		   double *Kzin, double *dzin, double *fb, double *ft,
		   int *ctp, int *cbm, int vcs, double *depth,
		   double *scale, double *w);
void do_vel_increment_3d(geometry_t *window);
void blend_vel(geometry_t *window, window_t *windat, win_priv_t *wincon, int mode, double *vel);
void calc_geostrophic(geometry_t *window, window_t *windat, win_priv_t *wincon);
void calc_geostrophic_c(geometry_t *window, window_t *windat, win_priv_t *wincon,
			int cs, double *vel, int mode);
void calc_geostrophic_obc(geometry_t *window, window_t *windat, win_priv_t *wincon,
			  open_bdrys_t *open, double *vel, int mode);
void rotate_vel(geometry_t *window, double *u, double *v);
void rerotate_vel(geometry_t *window, double *u, double *v);
void merge_thin_layers(int c, int zm1, double *vel, double *dz);
void vel_cen(geometry_t *window, window_t *windat, win_priv_t *wincon, double *u1,
	     double *u2, double *u, double *v, double *mag, double *dir, int mode);
void vel_grad(geometry_t *window, window_t *windat, win_priv_t *wincon,
	      double *u1, double *u2, double *ung, double *utg, int mode);
void tra_grad(geometry_t *window, window_t *windat, win_priv_t *wincon,
	      double *tr, double *ung, double *utg, int mode);
void vel_tan_3d(geometry_t *window, window_t *windat, win_priv_t *wincon);
double vel_c2e(geometry_t *geom, double *u, double *v, int e);

/*------------------------------------------------------------------*/
/* Turbulence closure routines                                      */
/*------------------------------------------------------------------*/
void closure_init(parameters_t *params, master_t *master);
void calcVzKz(geometry_t *window, window_t *windat, win_priv_t *wincon);
void bdry_closure(geometry_t *window, window_t *windat, win_priv_t *wincon);
void smooth_closure(geometry_t *window, window_t *windat, win_priv_t *wincon);
void closure_MY2(geometry_t *window, window_t *windat, win_priv_t *wincon);
void closure_constant(geometry_t *window, window_t *windat,
                      win_priv_t *wincon);
void closure_k_e_implicit(geometry_t *window, window_t *windat,
                          win_priv_t *wincon);
void closure_k_w_implicit(geometry_t *window, window_t *windat,
                          win_priv_t *wincon);
void k_e_1layer(geometry_t *window, window_t *windat, win_priv_t *wincon,
		double *wu, double *wv, int c, int cs);
void tridiagonal(geometry_t *window, int cs, int cb, double *C,
                 double *Cm1, double *Cp1, double *rhs, double *var,
                 double minvar);
void s_null(double aN, double aM, double *cmu, double *cmu_d);
void s_canutoA(double aN, double aM, double *cmu, double *cmu_d);
void s_canutoB(double aN, double aM, double *cmu, double *cmu_d);
void s_galperin(double aN, double aM, double *cmu, double *cmu_d);
void s_kantha_clayson(double aN, double aM, double *cmu, double *cmu_d);
void s_pom(double aN, double aM, double *cmu, double *cmu_d);
void s_munk_and(double aN, double aM, double *cmu, double *cmu_d);
void s_schum_gerz(double aN, double aM, double *cmu, double *cmu_d);
void s_eifler_schrim(double aN, double aM, double *cmu, double *cmu_d);

/*------------------------------------------------------------------*/
/* Dumpdata / dumpfile routines                                     */
/*------------------------------------------------------------------*/
dump_data_t *create_dumpdata(void);
dump_data_t *dumpdata_build(parameters_t *params, geometry_t *geom,
                            master_t *master, unsigned long ***flag,
                            double *layers, double **bathy);
void dumpdata_init(dump_data_t *dumpdata, geometry_t *geom,
                   master_t *master);
void dumpdata_init_geom(parameters_t *params, geometry_t *geom, dump_data_t *dumpdata);
void dumpdata_fill(geometry_t *geom, master_t *master,
                   dump_data_t *dumpdata);
void dumpdata_cleanup(dump_data_t *dumpdata, int mode);
int dump_open_us(parameters_t *params, char *name, int check);
void dump_close(int cdfid);
int dump_choose_by_time(parameters_t *params, int fid, double t);
int dump_choose_by_time_p(parameters_t *params, int fid, double t);
int dump_choose_by_time_m(master_t *master, int fid, double t);
int dump_choose_by_time_mom(master_t *master, int fid, double t);
int dump_choose_by_time_s(int fid, double t);
void read_grid_atts(parameters_t *params, int cdfid);
void read_mean_atts(master_t *master, int fid);
int get_nc_mode(dump_file_t *df);
int dump_re_read(master_t *master, 
		 int cdfid, int ti);
int dumpdata_read_us(geometry_t *geom, parameters_t *params, master_t *master,
		     dump_data_t *dumpdata, int cdfid, int ti);
int dumpdata_read_2d_us(dump_data_t *dumpdata, int id, char *name,
			double *p, int dump, int ni, int oset);
int dumpdata_read_3d_us(dump_data_t *dumpdata, int id, char *name,
			double *p, double **d, int dump, int ni, int ni3, int *map, int **rmap, int oset);
void create_df(dump_data_t *dumpdata, char *name, int mode);
int dump_init(sched_event_t *event);
void dumpfile_setup(master_t *master);
void dumpfile_resetup(master_t *master);
void dump_cleanup(sched_event_t *event, double t);
double dump_event(sched_event_t *event, double t);
void dump_snapshot(sched_event_t *event, double t);
void dump_f_snapshot(sched_event_t *event, int f, double t);
void dump_m_snapshot(master_t *master);
void dump_mp_snapshot(master_t *master);
void* dump_dispatch(void* data);
int dump_progress(sched_event_t *event);
void dump_windows_us(master_t *master, geometry_t **window, char *name, char *iname);
void *df_std_create(dump_data_t *dumpdata, dump_file_t *df);
void df_std_write(dump_data_t *dumpdata, dump_file_t *df, double t);
void df_std_close(dump_data_t *dumpdata, dump_file_t *df);
void *df_sp_create(dump_data_t *dumpdata, dump_file_t *df);
void df_sp_write(dump_data_t *dumpdata, dump_file_t *df, double t);
void df_sp_close(dump_data_t *dumpdata, dump_file_t *df);
void *df_ugrid_create(dump_data_t *dumpdata, dump_file_t *df);
void df_ugrid_write(dump_data_t *dumpdata, dump_file_t *df, double t);
void df_ugrid_close(dump_data_t *dumpdata, dump_file_t *df);
void df_ugrid_reset(dump_data_t *dumpdata, dump_file_t *df, double t);
void *df_ugrid3_create(dump_data_t *dumpdata, dump_file_t *df);
void df_ugrid3_write(dump_data_t *dumpdata, dump_file_t *df, double t);
int df_ugrid_get_varinfo(dump_data_t *dumpdata, dump_file_t *df,
			 char *name, df_ugrid_var_t *var);
void *df_parray_create(dump_data_t *dumpdata, dump_file_t *df);
void df_parray_write(dump_data_t *dumpdata, dump_file_t *df, double t);
void df_parray_close(dump_data_t *dumpdata, dump_file_t *df);
void *df_mom_create(dump_data_t *dumpdata, dump_file_t *df);
void df_mom_write(dump_data_t *dumpdata, dump_file_t *df, double t);
void df_mom_close(dump_data_t *dumpdata, dump_file_t *df);
int find_next_restart_record(dump_file_t *df, int cdfid, char *timevar, double tref);
int read3d(geometry_t *geom, int id, char *name, double *p, int dump,
            int nk, int nj, int ni);
/* Structure to describe each dump file time dep variable */
momgrid_t *convert_mom_grid(parameters_t *params, dump_data_t *dumpdata);
romsgrid_t *convert_roms_grid(parameters_t *params, dump_data_t *dumpdata);
void write_mom_grid(dump_data_t *dumpdata);
void write_roms_grid(dump_data_t *dumpdata);
void read_mom_grid(char *fname, int nx, int ny, double **x, double **y);
void read_roms_grid(char *fname, int nx, int ny, double **x, double **y, romsgrid_t *romsgrid);
void mom_grid_free(momgrid_t *momgrid, int mode);
void roms_grid_free(romsgrid_t *momgrid);
double *read_mom_layers(char *fname, int *nz);
double *read_bathy_from_mom(char *fname, int nx, int ny);		   
double *read_bathy_from_roms(char *fname, int nx, int ny, romsgrid_t *romsgrid);
void mask_bathy_from_roms(char *fname, parameters_t *params, romsgrid_t *romsgrid);

int hd_ts_multifile_get_index(int ntsfiles, timeseries_t **tsfiles,
                                cstring *names, char *var, int *varids);
timeseries_t *hd_ts_read(master_t *master, char *name, int check);
void hd_ts_free(master_t *master, timeseries_t *ts);
void hd_ts_check(master_t *master, timeseries_t *ts);
void timeseries_init(FILE * prmfd, master_t *master,geometry_t *geom,
		     geometry_t **window, dump_data_t *dumpdata);
void ts_resetup(master_t *master, FILE *fp);
void timeseries_init_w(master_t *master, geometry_t **window);
void timeseries_end(void);
timeseries_t **hd_ts_multifile_read(master_t *master, int nf,
                                    cstring * files);
int hd_ts_multifile_check(int ntsfiles, timeseries_t **tsfiles,
                          cstring * names, char *var, double tstart,
                          double tstop);
double hd_ts_multifile_eval_xy(int ntsfiles, timeseries_t **tsfiles,
                                int *varids, double t,
                                double x, double y);
double hd_ts_multifile_eval_xyz(int ntsfiles, timeseries_t **tsfiles,
                                int *varids, double t,
                                double x, double y, double z);
double hd_ts_multifile_eval_xy_by_name(int ntsfiles, timeseries_t **tsfiles,
                                cstring * names, char *var, double t,
                                double x, double y);
double hd_ts_multifile_eval_xyz_by_name(int ntsfiles, timeseries_t **tsfiles,
                                cstring * names, char *var, double t,
                                double x, double y, double z);
void hd_ts_multifile_eval_sparse(int ntsfiles, timeseries_t **tsfiles,
				 cstring * names, char *varname, double *var,
				 double t, int ns, int thio);
void hd_ts_multifile_eval_isparse(int ntsfiles, timeseries_t **tsfiles,
				  cstring * names, char *varname, double *var,
				  double t, int ns);
void hd_ts_multifile_eval(master_t *master, 
			  int ntsfiles, timeseries_t **tsfiles,
			  cstring * names, char *var, double *v, double t,
			  int *vec, int nvec, int mode);
void hd_vel_multifile_eval(master_t *master, 
			   int ntsfiles, timeseries_t **tsfiles,
			   cstring * names, char *var, double *v, double t,
			   int *vec, int nvec, int mode);
char *hd_ts_multifile_get_text_att(int ntsfiles, timeseries_t **tsfiles,
				   cstring * names, char *varname, char *attname);
void hd_trans_multifile_eval(master_t *master, 
			     int ntsfiles, timeseries_t **tsfiles,
			     cstring * names, char *var, double *v, double t,
			     int *vec, int nvec, int mode);
int hd_xyztoijk(void *data, double x, double y, double z, int *index,
                int *n, int *k);
int hd_grid_xytoij(master_t *master, double x, double y, int *i, int *j);
int hd_grid_xytoij_r(master_t *master, double x, double y, int *i, int *j);
int hd_grid_xyztoc(geometry_t *window, double x, double y, double z, int ci);
int hd_grid_xyztoc_w(geometry_t *window, double x, double y, double z);
int hd_grid_xyztoc_m(master_t *master, double x, double y, double z);
int hd_get_tracer_index_m(void *data, char *name);
int hd_get_tracer_index_w(void *data, char *name);

int hd_xyztoindex_m(void *data, double x, double y, double z,
                    int *cr, int *cs, int *cb);
int hd_xyztoindex_w(void *data, double x, double y, double z,
                    int *cr, int *cs, int *cb);
int xytoc(geometry_t *geom, xytoij_tree_t *xyij_tree, double x, double y, double *w1, double *w2);
long ztok(master_t *master, double z);
long ztoc_m(master_t *master, int c, double z);
long ztoc_w(geometry_t *window, int c, double z);
long ztoc(geometry_t *geom, int c, double z);
long ztoc_z(geometry_t *window, int c, double z, double *cz);
long ztoc_dz(geometry_t *window, int c, double z, double *dz, double *fz);
int ncu_dim_id(int cdfid, const char *name);
int ncu_var_id(int cdfid, const char *name);
void set_dims(int dimf, int *dims, int d1, int d2);
void set_count(int dimf, size_t *count, int d1, int d2);
landfillfn_t locate_landfill_function(const char *tag);
void dump_bathy_mask(dump_data_t *dumpdata, double bathyf);
void write_text_att(int cdfid, int varid, const char *name,
                    const char *text);
void dump_windows(master_t *master, geometry_t **window, char *name, char *iname);
void read_windows(geometry_t *geom, geometry_t **window, char *name);
void read_windows_us(geometry_t *geom, geometry_t **window, char *name);
void trans_check_dump(master_t *master, dump_data_t *dumpdata, char *trdata);
void hd_grid_interp_init(GRID_SPECS *gs, double *tr, char *method);
void hd_ts_grid_interp(master_t *master, timeseries_t *ts, char *varname,
		       double *tr, double t, int *vec, int nvec, char *method);
void hd_ts_grid_interp_multifile(master_t *master, timeseries_t **ts, int ntsfiles, int *varids,
				 double *tr, double t, int *vec, int nvec, char *method);
int find_closest_nonnan(dump_data_t *dumpdata, double **a, 
			int oi, int oj, int k, int *ci, int *cj);
int find_closest_w(geometry_t *geom, int oi, int oj, int *ci, int *cj);
int find_closest_b(parameters_t *params, double **bathy, double lat, 
		   double lon, int *ci, int *cj, int mode);
void cs_fill_land(dump_data_t *dumpdata);
double next_year(double time, char *unit);
double next_season(double time, char *unit, int *smon);
double prev_season(double time, char *unit, int *season, int *smon);
double next_month(double time, char *unit, int *smon);
double prev_month(double time, char *unit, int *smon);
double next_day(double time, char *unit, int *sday);
double prev_day(double time, char *unit, int *sday);
int yrday(int year, int mon, int day);
void df_std_reset(dump_data_t *dumpdata, dump_file_t *df, double t);
void df_simple_reset(dump_data_t *dumpdata, dump_file_t *df, double t);
void df_parray_reset(dump_data_t *dumpdata, dump_file_t *df, double t);
void df_sp_reset(dump_data_t *dumpdata, dump_file_t *df, double t);
void df_mom_reset(dump_data_t *dumpdata, dump_file_t *df, double t);
int ts_init(sched_event_t *event);
int get_glider_loc(master_t *master, ts_point_t *tslist, timeseries_t *loc_ts, double t);
void init_2way(master_t *master, geometry_t **window, window_t **windat, win_priv_t **wincon);

/*------------------------------------------------------------------*/
/* NetCDF convienience utilities                                    */
/*------------------------------------------------------------------*/
#if 0
void nc_s_writesub_1d(int fid, int varid,
                             size_t * start, size_t * count,
                             short *values);
#endif
void nc_l_writesub_3d(int fid, int varid,
		      size_t * start, size_t * count,
                             unsigned long ***v);
void nc_i_writesub_1d(int fid, int varid, size_t * start,
		      size_t * count, int *values);
void nc_i_writesub_2d(int fid, int varid, size_t * start,
		      size_t * count, int **values);
void nc_d_writesub_1d(int fid, int varid, size_t * start,
		      size_t * count, double *values);
void nc_d_writesub_2d(int fid, int varid, size_t * start,
		      size_t * count, double **values);
void nc_d_writesub_3d(int fid, int varid, size_t * start,
		      size_t * count, double ***values);
nc_type nc_get_default_type(int bpv);

/*------------------------------------------------------------------*/
/* Density routines                                                 */
/*------------------------------------------------------------------*/
double eos(double s, double t, double p);
void eos2(double s, double t, double p, double *dens, double *dens_0);
double lindensity_w(double s, double t, double p);
double quaddensity_w(double s, double t, double p);
void Set_lateral_BC_density_w(double *dens, int sgbpt, int *bpt, int *bin);
void density_w(geometry_t *window, window_t *windat, win_priv_t *wincon);
void density_m(master_t *master);
void compute_ref_density(master_t *master, geometry_t *window,
                         window_t *windat, win_priv_t *wincon);
void density_c(geometry_t *window, window_t *windat, win_priv_t *wincon,
	       int cs);
void density_gc(geometry_t *window, window_t *windat, win_priv_t *wincon,
		int cs, int *map);

/*------------------------------------------------------------------*/
/* Particle tracking routines                                       */
/*------------------------------------------------------------------*/
double ptrack_event(sched_event_t *event, double t);
int ptrack_init(sched_event_t *event);
void ptrack_cleanup(sched_event_t *event, double t);
void pt_params_init(master_t *master, FILE * fp);
void pt_update(master_t *master, geometry_t **window, window_t **windat, win_priv_t **wincon);
void pt_setup(master_t *master, geometry_t **window, window_t **windat, win_priv_t **wincon);
void ptrack_end(void);

/*------------------------------------------------------------------*/
/* Source / sink routines                                           */
/*------------------------------------------------------------------*/
void sourcesink_init(parameters_t *params, master_t *master,
                     geometry_t **window, window_t **windat,
                     win_priv_t **wincon);
void sourcesink(master_t *master);
void ss_tracer(geometry_t *window, window_t *windat, win_priv_t *wincon,
               int n, double *dtracer, double dt);
void ss_water(geometry_t *window, window_t *windat, win_priv_t *wincon);
void ss_momentum(geometry_t *window, window_t *windat, win_priv_t *wincon, int mode);
void sourcesink_reinit(master_t *master, geometry_t **window,
		       window_t **windat, win_priv_t **wincon);
void ref_depth_m(master_t *master, pss_t *pss, double *zlow, double *zhigh, int cs);

/*------------------------------------------------------------------*/
/* Relaxation routines                                              */
/*------------------------------------------------------------------*/
relax_info_t *relax_info_init(char *rname, char *tcname, double dt, int sz1, int sz2);
void relax_fill(geometry_t *window, relax_info_t *in, relax_info_t *out);

/*------------------------------------------------------------------*/
/* Standard boundary functions                                      */
/*------------------------------------------------------------------*/
void bf_hdstd_to_u1_init_m(master_t *master, open_bdrys_t *open,
                           bdry_details_t *d);
void bf_hdstd_to_u1_init_w(geometry_t *window, open_bdrys_t *open,
                           bdry_details_t *d, bdry_details_t *din);
void bf_hdstd_to_u1_m(geometry_t *geom, master_t *master,
                      open_bdrys_t *open, bdry_details_t *data);
double bf_hdstd_to_u1_w(geometry_t *window, window_t *windat,
                        win_priv_t *wincon, open_bdrys_t *open, double t,
                        int c, int cc, bdry_details_t *data);
void bf_hdstd_to_u1_t(master_t *master, open_bdrys_t *open_w,
                      bdry_details_t *data, geometry_t *window,
                      window_t *windat);
void bf_hdstd_to_u1av_m(geometry_t *geom, master_t *master,
                      open_bdrys_t *open, bdry_details_t *data);

void bf_hdstd_to_u2_init_m(master_t *master, open_bdrys_t *open,
                           bdry_details_t *d);
void bf_hdstd_to_u2_init_w(geometry_t *window, open_bdrys_t *open,
                           bdry_details_t *d, bdry_details_t *din);
void bf_hdstd_to_u2_m(geometry_t *geom, master_t *master,
                      open_bdrys_t *open, bdry_details_t *data);
double bf_hdstd_to_u2_w(geometry_t *window, window_t *windat,
                        win_priv_t *wincon, open_bdrys_t *open, double t,
                        int c, int cc, bdry_details_t *data);
void bf_hdstd_to_u2_t(master_t *master, open_bdrys_t *open_w,
                      bdry_details_t *data, geometry_t *window,
                      window_t *windat);
void bf_hdstd_to_u2av_m(geometry_t *geom, master_t *master,
                      open_bdrys_t *open, bdry_details_t *data);

void bf_uv_to_u1_init_m(master_t *master, open_bdrys_t *open,
                        bdry_details_t *d);
void bf_uv_to_u1_init_w(geometry_t *window, open_bdrys_t *open,
                        bdry_details_t *d, bdry_details_t *din);
void bf_uv_to_u1_m(geometry_t *geom, master_t *master, open_bdrys_t *open,
                   bdry_details_t *data);
double bf_uv_to_u1_w(geometry_t *window, window_t *windat,
                     win_priv_t *wincon, open_bdrys_t *open, double t,
                     int c, int cc, bdry_details_t *data);
void bf_uv_to_uav_init_w(geometry_t *window, open_bdrys_t *open,
			 bdry_details_t *data, bdry_details_t *data_in);
void bf_uv_to_u1av_m(geometry_t *geom, master_t *master, open_bdrys_t *open,
		     bdry_details_t *data);
double bf_uv_to_u1av_w(geometry_t *window, window_t *windat,
		       win_priv_t *wincon, open_bdrys_t *open, double t,
		       int c, int cc, bdry_details_t *data);
void bf_uv_to_u1_t(master_t *master, open_bdrys_t *open_w,
                   bdry_details_t *data, geometry_t *window,
                   window_t *windat);

void bf_uv_to_u2_init_m(master_t *master, open_bdrys_t *open,
                        bdry_details_t *d);
void bf_uv_to_u2_init_w(geometry_t *window, open_bdrys_t *open,
                        bdry_details_t *d, bdry_details_t *din);

void bf_uv_to_u2_m(geometry_t *geom, master_t *master,
                   open_bdrys_t *open, bdry_details_t *data);
double bf_uv_to_u2_w(geometry_t *window, window_t *windat,
                     win_priv_t *wincon, open_bdrys_t *open, double t,
                     int c, int cc, bdry_details_t *data);
void bf_uv_to_u2av_m(geometry_t *geom, master_t *master,
		     open_bdrys_t *open, bdry_details_t *data);
double bf_uv_to_u2av_w(geometry_t *window, window_t *windat,
		       win_priv_t *wincon, open_bdrys_t *open, double t,
                     int c, int cc, bdry_details_t *data);
void bf_uv_to_u2_t(master_t *master, open_bdrys_t *open_w,
                   bdry_details_t *data, geometry_t *window,
                   window_t *windat);

void bf_uvav_to_u1av_init_m(master_t *master, open_bdrys_t *open,
			    bdry_details_t *d);
void bf_uvav_to_u1av_init_w(geometry_t *window, open_bdrys_t *open,
			    bdry_details_t *d, bdry_details_t *din);
void bf_uvav_to_u1av_m(geometry_t *geom, master_t *master, 
		       open_bdrys_t *open, bdry_details_t *data);
double bf_uvav_to_u1av_w(geometry_t *window, window_t *windat,
			 win_priv_t *wincon, open_bdrys_t *open, double t,
			 int c, int cc, bdry_details_t *data);
void bf_uvav_to_u1av_t(master_t *master, open_bdrys_t *open_w,
		       bdry_details_t *data, geometry_t *window,
		       window_t *windat);

void bf_uvav_to_u2av_init_m(master_t *master, open_bdrys_t *open,
			    bdry_details_t *d);
void bf_uvav_to_u2av_init_w(geometry_t *window, open_bdrys_t *open,
			    bdry_details_t *d, bdry_details_t *din);
void bf_uvav_to_u2av_m(geometry_t *geom, master_t *master, 
		       open_bdrys_t *open, bdry_details_t *data);
double bf_uvav_to_u2av_w(geometry_t *window, window_t *windat,
			 win_priv_t *wincon, open_bdrys_t *open, double t,
			 int c, int cc, bdry_details_t *data);
void bf_uvav_to_u2av_t(master_t *master, open_bdrys_t *open_w,
		       bdry_details_t *data, geometry_t *window,
		       window_t *windat);

void bf_u1_flow_init_m(master_t *master, open_bdrys_t *open,
                       bdry_details_t *d);
void bf_u1_flow_init_w(geometry_t *window, open_bdrys_t *open,
                       bdry_details_t *d, bdry_details_t *din);
void bf_flow_m(geometry_t *geom, master_t *master, open_bdrys_t *open,
               bdry_details_t *d);
double bf_u1_flow_w(geometry_t *window, window_t *windat,
                    win_priv_t *wincon, open_bdrys_t *open, double t,
                    int c, int cc, bdry_details_t *data);
void bf_flow_t(master_t *master, open_bdrys_t *open_w,
               bdry_details_t *data, geometry_t *window, window_t *windat);

void bf_u2_flow_init_m(master_t *master, open_bdrys_t *open,
                       bdry_details_t *d);
void bf_u2_flow_init_w(geometry_t *window, open_bdrys_t *open,
                       bdry_details_t *d, bdry_details_t *din);
double bf_u2_flow_w(geometry_t *window, window_t *windat,
                    win_priv_t *wincon, open_bdrys_t *open, double t,
                    int c, int cc, bdry_details_t *data);

void bf_use_eqn_init_m(master_t       *master, 
		       open_bdrys_t   *open, 
		       bdry_details_t *data);
void bf_use_eqn_m(geometry_t *geom, master_t *master,
		  open_bdrys_t *open, bdry_details_t *data);
void bf_use_eqn_trans(master_t *master, open_bdrys_t *open_w,
		      bdry_details_t *data, geometry_t *window,
		      window_t *windat);
double bf_use_eqn_w(geometry_t *window, window_t *windat,
		    win_priv_t *wincon, open_bdrys_t *open, double t,
		    int c, int cc, bdry_details_t *data);

void bf_ts_free(master_t *master, bdry_details_t *data);
void bf_c2cc_free(geometry_t *window, bdry_details_t *data);
void bf_void_free(geometry_t *window, bdry_details_t *data);
void bf_flow_free(master_t *master, bdry_details_t *data);
void bf_eqn_free(master_t *master, bdry_details_t *data);

/*------------------------------------------------------------------*/
/* Miscellaneous                                                    */
/*------------------------------------------------------------------*/
void hd_silent_warn(const char *s, ...);
void hd_warn(const char *s, ...);
void hd_error(const char *s, ...);
void hd_quit(const char *s, ...);
void hd_quit_and_dump(const char *s, ...);
void step_end_2d_init(master_t *master);
void step_end_2d(master_t *master);
void step_end_3d_init(master_t *master);
void step_end_3d(master_t *master);
void custdata_init(master_t *master);
void custdata_end(master_t *master);
void process_args(int argc, char *argv[]);
void readname(char *s, char *prompt);
double shuman(geometry_t *window, double *a, int c);
void shuman_3d(geometry_t *window, double *var, double v);
double cvol1(master_t *master, double *a, int c, int edge);
void smooth_w(geometry_t *window, double *A, double *B,
	      int *ctp, int nctp, int n, double lim, int edgef);
void shuman_smooth(geometry_t *window, double *a, int cl);
void shapiro(geometry_t *window, double *a, double *buf, int *vec, 
	     int nvec, int order, int mode, int dir);
void shapiro_smooth(geometry_t *window, double *a, int cl, int mode);
void smooth3(master_t *master, double *A, int *ctp, int nctp, int sz, int edge);
void smooth3e(master_t *master, double *A, int *ctp, int nctp, int sz);
void smooth(master_t *master, double *A, int *ctp, int nctp);
void free1ds(short *p);
int ANY(int var, int array[], int ns);
int ANY0(int var, int array[], int ns);
char *fv_get_filename(char *ofn, char *buf);
char *fv_get_varname(char *oprm, char *tag, char *buf);
int is_true(const char *tag);
void printflags(FILE * out, unsigned long f);
int find_non_nan(double **a, int nx, int ny, int oi, int oj, int *ci, int *cj);
/*UR-ADDED */
void free_stack_add(void* object,void (*free) (void*) );
 custom_function_t* custom_stack_remove(char* fnc);
char* custom_stack_add(custom_function_t* fnc);
 void custom_init(hd_data_t* hdata);
void custom_step(hd_data_t* hdata);
double intp(double a, double b, int xs, int xe, int x);
double intpf(double a, double b, double xs, double xe, double x);

/*------------------------------------------------------------------*/
/* ginterface                                                       */
/*------------------------------------------------------------------*/
void i_set_error(void* hmodel, int col, int errorf, char *text);
int i_get_error(void* hmodel, int col);
void ginterface_moonvars(void *hmodel, int c,
			 double *mlon, double *mlat,
			 double *dist_earth_sun, double *dist_moon_earth,
			 double *lunar_angle, double *sun_angle,double *moon_phase,double *lunar_dec);
double ginterface_get_cloud(void *hmodel, int c);


/*------------------------------------------------------------------*/
/* Ecology                                                          */
/*------------------------------------------------------------------*/
#if defined(HAVE_ECOLOGY_MODULE)
void eco_step(geometry_t *window);
void eco_set_tracer_defaults(tracer_info_t *tracer, char *trname,
			     char *defname, ecology *e);
void eco_read_tr_atts(tracer_info_t *tr, FILE *fp, char *keyname);
void read_ecology(parameters_t *params, FILE *fp, int *ntr);
int ecology_autotracer_2d(FILE *fp, int do_eco, char *eco_vars, char *eco_defs,
			  void *e, tracer_info_t *trinfo, int ntr, int tn);
int ecology_autotracer_3d(FILE *fp, int do_eco, char *eco_vars, char *eco_defs,
			  void *e, tracer_info_t *trinfo, int ntr, int tn);
int ecology_autotracer_sed(FILE *fp, int do_eco, char *eco_vars, char *eco_defs,
			   void *e, tracer_info_t *trinfo, int ntr, int tn);
int count_eco_3d(char *eco_vars);
int count_eco_2d(char *eco_vars);
void eco_write_tr_atts(tracer_info_t *tr, FILE *fp, int n);
int ecology_autotracer_write(master_t *master, FILE *op, int tn);
int is_eco_var(char *trname);
int eco_get_obc(tracer_info_t *tr);
int get_eco_var_type(char *trname, char *eco_defs);
void trans_write_eco(char *eco_vars, char *eco_defs, double ecodt, ecology *e, FILE *fp);
int eco_get_nstep(ecology *e);
void *ecology_pre_build(char *eco_vars, char *eco_defs, FILE *prmfd);
void eco_set_omp_num_threads(ecology *e, int n);
double get_parameter_value(ecology* e, char* s);
double try_parameter_value(ecology* e, char* s);
#endif

/*------------------------------------------------------------------*/
/* Waves                                                            */
/*------------------------------------------------------------------*/
#if defined(HAVE_WAVE_MODULE)
void wave_interface_step(geometry_t *window);
#endif
double get_bv_mono(geometry_t *window, window_t *windat, win_priv_t *wincon,
		   int c, int cs);

/*------------------------------------------------------------------*/
/* Sediments                                                        */
/*------------------------------------------------------------------*/
#if defined(HAVE_SEDIMENT_MODULE)
sediment_t *sed_init(FILE * prmfd, void *model);
void hd2sed(sediment_t *sediment, sed_column_t *sm, int c);
void sed2hd(sediment_t *sediment, sed_column_t *sm, int c);
void sed_step(geometry_t *window);
void sed_cleanup(sediment_t *sediment);

void sed_set_grid(geometry_t *geom, geometry_t **window);
void eco_sed_init(FILE * fp, geometry_t *window);
void read_sediments(parameters_t *params, FILE *fp, int *ntr);
void sed_set_tracer_defaults(tracer_info_t *tracer, char *trname, char *defname);
void sed_read_tr_atts(tracer_info_t *tr, FILE *fp, char *keyname);
void sed_set_att(tracer_info_t *tr, char *key, void *value);
void sed_tracer_custom(master_t *master);
int count_sed_classes(char *sed_vars);
int count_sed_3d();
int count_sed_2d();
int sediment_autotracer_3d(FILE *fp, int do_sed, char *sed_vars, char *sed_defs, 
			    tracer_info_t *trinfo, int ntr, int tn);
int sediment_autotracer_2d(FILE *fp, int do_sed, char *sed_vars, char *sed_defs, 
			    tracer_info_t *trinfo, int ntr, int tn);
int sediment_autotracer_sed(FILE *fp, int do_sed, char *sed_vars, char *sed_defs, 
			     tracer_info_t *trinfo, int ntr, int tn);
int sediment_autotracer_write(master_t *master, FILE *op, int tn);
void sed_write_tr_atts(tracer_info_t *tr, FILE *fp, int n);
void print_tr_sed_atts(tracer_info_t *tr);
int sed_get_obc(tracer_info_t *tracer);
int is_sed_var(char *trname);
void trans_write_sed(parameters_t *params, sediment_t *sediment, FILE *fp);
#endif

/*------------------------------------------------------------------*/
/* Tracerstats                                                      */
/*------------------------------------------------------------------*/
#if defined(HAVE_TRACERSTATS_MODULE)
void tracerstats_step(geometry_t *window);
void tracerstats_prestep(geometry_t *window, int step);
#endif

#endif
